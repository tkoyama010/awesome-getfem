
#!/bin/sh
# update transifex pot and local po files
set -ex
LOCAL_PYTHON_PATH="/home/runner/.local/bin"
# pull po files from transifex
cd `dirname $0`
$LOCAL_PYTHON_PATH/sphinx-intl create-transifexrc
$LOCAL_PYTHON_PATH/sphinx-build -T -b gettext ../doc pot
$LOCAL_PYTHON_PATH/sphinx-intl update-txconfig-resources -p pot -d .
cat .tx/config
$LOCAL_PYTHON_PATH/tx push -s --skip
rm -Rf ja pt_BR
$LOCAL_PYTHON_PATH/tx pull -l ja,pt_BR
git checkout .tx/config
