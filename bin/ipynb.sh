#!/bin/bash

export PYTHONPATH=$PYTHONPATH:$HOME/bin

pushd /home/oceano/notebooks
  ipython3 notebook --matplotlib=inline
popd
