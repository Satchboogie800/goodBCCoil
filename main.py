"""
High-pass birdcage coil (16 legs) for Ansys HFSS via PyAEDT
============================================================

Builds:
  * zero-thickness PEC cylinder sheet, cut into 16 legs + 2 end rings
  * 16 capacitor gaps per end ring (caps sit between the legs)
  * 2 feed positions in the bottom ring (default 0 deg and 90 deg) using the
    split-capacitor match:   ring --|C1|-- island --|C2 || Port|-- ring
    (C2 and the port share the C2 gap, each over half the ring width)
  * optional PEC shield, air region with radiation boundary, setup + sweep

Two modes (CAPS_AS_PORTS):
  False -> every capacitor is a lumped RLC (final, "physical" model)
  True  -> every capacitor is a lumped port (for Circuit co-simulation /
           tuning; afterwards write the tuned values back and switch to False)

Written for ansys-aedt-core (PyAEDT >= 0.11, `from ansys.aedt.core import Hfss`).
NOT tested against a running AEDT here - argument names have changed between
PyAEDT releases, so if a call fails, check its signature in the API docs for
your installed version.

Geometry convention
-------------------
Coil axis = z, coil from z = 0 to z = L_COIL.
Leg-gap boxes are centred at k*22.5 deg  -> legs sit at 11.25 + k*22.5 deg,
end-ring capacitors sit at k*22.5 deg (between the legs).
All lengths in mm.
"""

import math
from ansys.aedt.core import Hfss

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
N_LEGS = 16
R = 150.0           # coil radius
L_COIL = 250.0      # total coil length incl. end rings
W_RING = 20.0       # end-ring width (axial)
W_LEG = 20.0        # leg width (arc length)
G_CAP = 2.5         # capacitor gap length (along current direction)
S_ISLAND = 4.0      # conductor island between C1 and C2 at the feeds
FEED_ANGLES = [0.0, 90.0]   # feed positions in the bottom ring (deg)

USE_SHIELD = True
R_SHIELD = 180.0
L_SHIELD = 350.0

F0_MHZ = 127.7       # 3 T proton frequency; 297.2 for 7 T
C_RING_PF = 20.0     # start value - replace with your tuned value
CAPS_AS_PORTS = False

AEDT_VERSION = "2025.1"
NON_GRAPHICAL = False

# ----------------------------------------------------------------------------
# Derived quantities
# ----------------------------------------------------------------------------
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
    project="Birdcage",
    design="HP_Birdcage_16",
    solution_type="Modal",
    version=AEDT_VERSION,
    non_graphical=NON_GRAPHICAL,
    new_desktop=True,
)
hfss.modeler.model_units = "mm"

# Capacitor values as design variables -> tunable with Optimetrics
hfss["C_ring"] = f"{C_RING_PF}pF"
hfss["C1"] = "1.2*C_ring"   # series C1*C2/(C1+C2) = C_ring at the start
hfss["C2"] = "6*C_ring"     # larger C2 -> lower impedance seen by the port


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def cylinder_sheet(name, r, z0, h):
    """Open cylindrical sheet (no caps) of radius r from z0 to z0+h."""
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
    if CAPS_AS_PORTS:
        hfss.lumped_port(
            assignment=sheet.name,
            create_port_sheet=False,
            integration_line=line,
            impedance=50,
            name=f"P_{sheet.name}",
            renormalize=True,
        )
        return
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
# Shield, air region, radiation boundary
# ----------------------------------------------------------------------------
if USE_SHIELD:
    shield = cylinder_sheet("Shield", R_SHIELD, (L_COIL - L_SHIELD) / 2, L_SHIELD)
    hfss.assign_perfecte_to_sheets(shield.name, name="PEC_Shield")

# Head model: import / place it here (e.g. hfss.modeler.import_3d_cad(...))
# and assign tissue materials before solving.

pad = 300  # mm; keep the boundary roughly lambda/4 away from the coil
region = hfss.modeler.create_air_region(
    x_pos=pad, y_pos=pad, z_pos=pad, x_neg=pad, y_neg=pad, z_neg=pad,
    is_percentage=False,
)
hfss.assign_radiation_boundary_to_objects(region.name, name="Rad")

# ----------------------------------------------------------------------------
# Mesh refinement on gaps
# ----------------------------------------------------------------------------
hfss.mesh.assign_length_mesh(
    [s.name for s in gap_sheets], maximum_length=G_CAP / 2, name="GapMesh"
)

# ----------------------------------------------------------------------------
# Setup + sweep
# ----------------------------------------------------------------------------
setup = hfss.create_setup("Setup1")
setup.props["Frequency"] = f"{F0_MHZ}MHz"
setup.props["MaximumPasses"] = 15
setup.props["MaxDeltaS"] = 0.01
setup.update()
setup.create_frequency_sweep(
    unit="MHz",
    name="Sweep",
    start_frequency=0.8 * F0_MHZ,
    stop_frequency=1.2 * F0_MHZ,
    sweep_type="Interpolating",
)

# Quadrature drive (ideal hybrid). Flip the sign of 90deg if B1+ comes out
# weak in the centre (wrong rotation sense for your B0 direction).
# Signature of edit_sources differs between versions - check yours:
# hfss.edit_sources({"Feed1:1": ("1W", "0deg"), "Feed2:1": ("1W", "90deg")})

hfss.save_project()
hfss.release_desktop(close_projects=False, close_desktop=False)