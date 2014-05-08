#!/bin/bash

export PYTHONPATH=$PYTHONPATH:$HOME/bin

pushd /home/oceano/notebooks
  ipython notebook
popd
