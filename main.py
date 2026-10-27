import math
from ansys.aedt.core import Hfss
# Initial Parameters
N_LEGS = 16
R = 150.0           # coil radius
L_COIL = 250.0      # total coil length incl. end rings
W_RING = 20.0       # end-ring width (axial)
W_LEG = 20.0        # leg width (arc length)
G_CAP = 4         # capacitor gap length (along current direction)
S_ISLAND = 4.0      # conductor island between C1 and C2 at the feeds
FEED_ANGLES = [0.0, 90.0]   # feed positions in the bottom ring (deg)

F0_MHZ = 64       # (~T*42MHz)
C_RING_PF = 105     # Initial tuning capacitor value in pF

AEDT_VERSION = "2023.2"


# Derived Parameters
DPHI = 360.0 / N_LEGS
# width of the leg-gap box (chord), see derivation: b = 2R sin(pi/N - w/2R)
B_BOX = 2 * R * math.sin(math.pi / N_LEGS - W_LEG / (2 * R))
L_LEGS = L_COIL - 2 * W_RING

assert S_ISLAND + 2 * G_CAP < B_BOX, "feed region wider than the leg gap"

def pt(theta_deg, y, z):
    """Point on the cylinder at frame angle theta, offset y (chord dir.), height z."""
    phi = math.radians(theta_deg) + math.asin(y / R)
    return [R * math.cos(phi), R * math.sin(phi), z]


# ----------------------------------------------------------------------------
# Launch HFSS
# ----------------------------------------------------------------------------
hfss = Hfss(
    project="Birdcage_2",
    design="HP_Birdcage_16",
    solution_type="Modal",
    version=AEDT_VERSION,
    non_graphical=False,
    new_desktop=False,
)
hfss.modeler.model_units = "mm"

# Capacitor values as design variables -> tunable with Optimetrics
hfss["C_ring"] = f"{C_RING_PF}pF"
hfss["n"] = "0.3" # ratio of C1 and C2, use this for matching the ports impedance
hfss["C1"] = "1.2*C_ring"   # C1 and C2 are adjusted by n (relation changes the port impedance,
hfss["C2"] = "C1*n"     #      the total value of C1 and C2 stays equal to C_ring)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def cylinder_sheet(name, r, z0, h):
    circ = hfss.modeler.create_circle(
        "XY", [0, 0, z0], r, is_covered=False, name=name
    )
    hfss.modeler.sweep_along_vector(circ, [0, 0, h])
    return circ


def half_box(name, theta, y1, y2, z0, h):
    """Box cutting the cylinder only on one side (x > 0 in the local frame)."""
    box = hfss.modeler.create_box([0, y1, z0], [1.5 * R, y2 - y1, h], name=name)
    hfss.modeler.rotate(box, "Z", theta)
    return box


def gap_sheet(name, theta, y1, y2, z0, h):
    """Curved sheet exactly filling a gap between local offsets y1..y2."""
    pts = [pt(theta, y, z0) for y in (y1, 0.5 * (y1 + y2), y2)]
    arc = hfss.modeler.create_polyline(pts, segment_type="Arc", name=name)
    hfss.modeler.sweep_along_vector(arc, [0, 0, h])
    return arc


def int_line(theta, y1, y2, z):
    """Integration line across the gap (along the ring current)."""
    return [pt(theta, y1, z), pt(theta, y2, z)]


def assign_cap(sheet, line, cap_var):
    bnd = hfss.assign_lumped_rlc_to_sheet(
        sheet.name,
        start_direction=line,
        rlc_type="Parallel",
        capacitance=C_RING_PF * 1e-12,
        name=f"RLC_{sheet.name}",
    )
    # link to the design variable so Optimetrics can tune it
    bnd.props["Capacitance"] = cap_var
    bnd.update()


# ----------------------------------------------------------------------------
# Conductor: cylinder minus leg gaps minus capacitor gaps
# ----------------------------------------------------------------------------
coil = cylinder_sheet("Coil", R, 0, L_COIL)
tools = []

# 8 boxes through the axis -> 16 leg gaps (box rotated in 360/16 steps!)
for k in range(N_LEGS // 2):
    box = hfss.modeler.create_box(
        [-1.5 * R, -B_BOX / 2, W_RING], [3 * R, B_BOX, L_LEGS], name=f"LegCut{k}"
    )
    hfss.modeler.rotate(box, "Z", k * DPHI)
    tools.append(box)

gap_sheets = []          # for mesh refinement
cap_jobs = []            # (sheet, integration line, capacitor variable)
port_jobs = []           # (sheet, integration line, port name)

rings = {"Bot": 0.0, "Top": L_COIL - W_RING}
for ring, z0 in rings.items():
    zmid = z0 + W_RING / 2
    for k in range(N_LEGS):
        th = k * DPHI
        is_feed = ring == "Bot" and any(abs(th - f) < 1e-6 for f in FEED_ANGLES)

        if not is_feed:
            y1, y2 = -G_CAP / 2, G_CAP / 2
            tools.append(half_box(f"Cut_{ring}{k}", th, y1, y2, z0 - 1, W_RING + 2))
            s = gap_sheet(f"C_{ring}{k}", th, y1, y2, z0, W_RING)
            cap_jobs.append((s, int_line(th, y1, y2, zmid), "C_ring"))
            gap_sheets.append(s)
            continue

        # ---- feed: C1 gap, island, C2 gap (C2 and port side by side) ----
        a1, a2 = -(S_ISLAND / 2 + G_CAP), -S_ISLAND / 2      # C1 gap
        b1, b2 = S_ISLAND / 2, S_ISLAND / 2 + G_CAP          # C2 gap
        tools.append(half_box(f"Cut_C1_{k}", th, a1, a2, z0 - 1, W_RING + 2))
        tools.append(half_box(f"Cut_C2_{k}", th, b1, b2, z0 - 1, W_RING + 2))

        s1 = gap_sheet(f"C1_{k}", th, a1, a2, z0, W_RING)
        cap_jobs.append((s1, int_line(th, a1, a2, zmid), "C1"))

        # port on the outer half of the ring, C2 on the inner half
        sp = gap_sheet(f"Feed_{k}", th, b1, b2, z0, W_RING / 2)
        port_jobs.append((sp, int_line(th, b1, b2, z0 + W_RING / 4), f"Feed{len(port_jobs) + 1}"))

        s2 = gap_sheet(f"C2_{k}", th, b1, b2, z0 + W_RING / 2, W_RING / 2)
        cap_jobs.append((s2, int_line(th, b1, b2, z0 + 3 * W_RING / 4), "C2"))

        gap_sheets += [s1, sp, s2]

hfss.modeler.subtract(coil, tools, keep_originals=False)
hfss.assign_perfecte_to_sheets(coil.name, name="PEC_Coil")

# ----------------------------------------------------------------------------
# Capacitors and ports
# ----------------------------------------------------------------------------
for sheet, line, var in cap_jobs:
    assign_cap(sheet, line, var)

for sheet, line, pname in port_jobs:
    hfss.lumped_port(
        assignment=sheet.name,
        create_port_sheet=False,
        integration_line=line,
        impedance=50,
        name=pname,
        renormalize=True,
    )

# ----------------------------------------------------------------------------
# Mesh refinement on gaps
# ----------------------------------------------------------------------------
hfss.mesh.assign_length_mesh(
    [s.name for s in gap_sheets], maximum_length=G_CAP / 2, name="GapMesh"
)

hfss.save_project()
hfss.release_desktop(close_projects=False, close_desktop=False)