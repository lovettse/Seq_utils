PATH="$PATH:/data/slovett/scripts_organized/working"

MYFILES=`ls -F /data/slovett/scripts_organized/working | grep "\/" `
for FOLDER in $MYFILES; do
  PATH="$PATH:/data/slovett/scripts_organized/working/$FOLDER"
done
export PATH
