#!/bin/bash

# Copy files over to cluster

# Data storage path
OLD_PATH_D="/Users/jsquire/Documents/MRIDSS/MRIDSS/Data/"
NEW_PATH_D="/p/mridss/Data/"
# Base directory path
OLD_PATH="/Users/jsquire/Documents/MRIDSS/MRIDSS/"
NEW_PATH="/u/jsquire/MRIDSS/MRIDSS/"

sed -i "" "s;$OLD_PATH_D;$NEW_PATH_D;g" MRIDSS/General_Definitions.h
sed -i "" "s;$OLD_PATH;$NEW_PATH;g" MRIDSS/General_Definitions.h
echo "Changed path ${OLD_PATH} to ${NEW_PATH}"

# Only copy over source files
scp MRIDSS/main.cc portal:~/MRIDSS/MRIDSS/
scp MRIDSS/General_Definitions.h portal:~/MRIDSS/MRIDSS/
scp -r MRIDSS/Models portal:~/MRIDSS/MRIDSS/
scp -r MRIDSS/Auxiliary portal:~/MRIDSS/MRIDSS/
scp -r MRIDSS/Integrators portal:~/MRIDSS/MRIDSS/
# Matlab routines
# scp -r Matlab portal:~/MRIDSS/

echo ""
echo "Files updated!"
echo ""

sed -i "" "s;$NEW_PATH;$OLD_PATH;g" MRIDSS/General_Definitions.h
sed -i "" "s;$NEW_PATH_D;$OLD_PATH_D;g" MRIDSS/General_Definitions.h

echo "Changed path ${NEW_PATH} to ${OLD_PATH}"

