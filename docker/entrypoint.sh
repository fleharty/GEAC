#!/bin/bash
# Dispatch: "cohort" launches the cohort Explorer; everything else goes to geac.
if [ "$1" = "cohort" ]; then
    shift
    exec streamlit run /app/geac_explorer.py "$@"
else
    exec geac "$@"
fi
