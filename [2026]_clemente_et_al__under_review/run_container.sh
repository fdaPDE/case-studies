#!/bin/bash

nonlinear=$(pwd)
docker run --rm -v $nonlinear:/root/ -ti aldoclemente/fdapde-docker /bin/bash
