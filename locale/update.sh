
#!/bin/sh
# update transifex pot and local po files
set -ex
# pull po files from transifex
cd `dirname $0`
sphinx-intl create-transifexrc
sphinx-build -T -b gettext ../ pot
sphinx-intl update-txconfig-resources -p pot -d .
cat .tx/config
tx push -s --skip
rm -Rf ja
tx pull -l ja
git checkout .tx/config
