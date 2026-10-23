"""Sample definitions for the H->ZZ*->4l analysis (ATLAS Open Data 2025 education release)."""

RELEASE = "2025e-13tev-beta"
SKIM = "exactly4lep"
LUMI_FB = 36.6

# Grouping follows the official ATLAS Open Data HZZ notebook.
SAMPLES = {
    "data": ["data"],
    "higgs": [345060, 346228, 346310, 346311, 346312, 346340, 346341, 346342],
    "zz": [700600],
    "other": [410470, 410155, 410218, 410219, 412043, 364243, 364242, 364246,
              364248, 700320, 700321, 700322, 700323, 700324, 700325],
}
