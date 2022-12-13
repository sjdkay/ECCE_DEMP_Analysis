#! /bin/bash

### Stephen JD Kay, University of Regina
### 19/10/22 - The script simply rebuilds all of the analysis modules 5on100, 5on41 and 10on100
### This scrip should be run any time any of the three modules is changed

WorkDir="/home/fun4all/work" # This should be the only path you actually need to adjust
CodeDir5on100="${WorkDir}/ECCE_DEMP_Analysis/ECCE_DEMP_Ana"
BuildDir5on100="${WorkDir}/ECCE_DEMP/build"
CodeDir5on41="${WorkDir}/ECCE_DEMP_Analysis/ECCE_DEMP5on41_Ana"
BuildDir5on41="${WorkDir}/ECCE_DEMP5on41/build"
CodeDir10on100="${WorkDir}/ECCE_DEMP_Analysis/ECCE_DEMP10on100_Ana"
BuildDir10on100="${WorkDir}/ECCE_DEMP10on100/build"

if [ ! -d "${BuildDir5on100}" ] || [ ! -d "${BuildDir5on41}" ] || [ ! -d "${BuildDir10on100}" ]; then
    echo "This script is only intended to rebuild the re-analysis modules, please follow the instructions in the README for each module and install it first."
    exit 1
fi

cd "${BuildDir5on100}"
cp "${CodeDir5on100}/ECCE_DEMP."* "${BuildDir5on100}/"
make install
echo""
echo "5on100 done, rebuilding 5on41 next"
echo""

cd "${BuildDir5on41}"
cp "${CodeDir5on41}/ECCE_DEMP5on41."* "${BuildDir5on41}/"
make install
echo""
echo "5on41 done, rebuilding 10 on 100 next"
echo""

cd "${BuildDir10on100}"
cp "${CodeDir10on100}/ECCE_DEMP10on100."* "${BuildDir10on100}/"
make install
echo""
echo"Done rebuilding 10on100, all done! Check for any errors!"
echo""
exit 0


