#!/bin/bash

root -l -b -q makeAllAccEffMaps.C+'(1)'
cp -r FilesAccxEff_pt_SizeDoubled_centBins_NoWeights FilesAccxEff_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights
root -l -b -q makeAllAccEffMaps.C+'(0)'
