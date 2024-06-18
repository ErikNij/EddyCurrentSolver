#!/bin/bash

wArgs='-Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -Wfatal-errors'
oArgs='-m64 -O3'

g++ $oArgs $wArgs -o cppVirtual cppVirtual.C
