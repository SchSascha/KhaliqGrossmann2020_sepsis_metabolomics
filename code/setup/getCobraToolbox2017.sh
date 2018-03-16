#!/bin/bash
# This cript fetches the current version of the COBRA toolbox and places it in the directory "../cobratoolbox/". Also initializes submodules.
cd ..
git clone https://github.com/opencobra/cobratoolbox.git cobratoolbox
cd cobratoolbox
git submodule update --init #prevents fail of initCobraToolbox at getDistributedModelFolder()
