#!/bin/bash
# Dispatch: "explorer" launches the Streamlit app; everything else goes to geac.
if [ "$1" = "explorer" ]; then
    shift
    exec streamlit run /app/geac_explorer.py "$@"
else
    exec geac "$@"
fi
