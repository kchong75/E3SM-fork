here=`pwd`

# Contents of this scripts are based on structions at https://gist.github.com/slint/2263b2212743a68b2851b2cb7865675e
#---------------------------------------------------------------------------------------------------
# Step 1. Create an access token
#---------------------------------
# Create a token at "https://zenodo.org/account/settings/applications/tokens/new/"
# with the "deposit:write" and "deposit:actions" scopes. Once the token is created, save it as
# the variable TOKEN below
 
TOKEN="1234567"

#---------------------------------------------------------------------------------------------------
# Step 2. Start a new upload via the web interface
#--------------------------------------------------
# At the Zenodo website, start a new upload, fill in the minimum metadata (title, authors,
# description, access rights and license) and click "Save".

#---------------------------------------------------------------------------------------------------
# Step 3. Get the bucket number of the new upload
#--------------------------------------------------
# After the new upload has been started and the basic info saved, 
# look into your browser's URL box, find the deposit ID that comes after "https://zenodo/deposit/".
# Copy and save the deposit ID to the variable DEPOSIT_ID below, then 
# run the curl command with the output redirected to a txt file.

#DEPOSIT_ID="10407375"
#curl "https://zenodo.org/api/deposit/depositions/${DEPOSIT_ID}?access_token=${TOKEN}" >curl.out.log
#exit

#--

# in curl.out.log, look for the word "bucket". You will find something like
#
# "links": { "bucket": "https://zenodo.org/api/files/568377dd-daf8-4235-85e1-a56011ad454b", ... }
#
# copy the long string that comes after "files/" then save it to the variable BUCKET below.

BUCKET="08e5f45f-6efc-4783-9736-5a0e0559e0c5"

#---------------------------------------------------------------------------------------------------
# Step 4. Upload files
#--------------------------------------------------
# Note: to use the curl command below, we need the variables TOKEN and BUCKET declared above (but not DEPOSIT_ID).
#
# The curl command below upload the local file ${FILEPATH}/${FILENAME} as ${FILENAME} on Zenodo

FILEPATH="/pscratch/sd/h/huiwan/cflx/v1_papers/for_upload/"

cd ${FILEPATH}

nfile=0
for FILENAME in `ls  *.tar.gz* ` ;do

    echo =========
    echo Uploading $FILENAME
    echo
    curl --upload-file "${FILEPATH}/${FILENAME}" "https://zenodo.org/api/files/${BUCKET}/${FILENAME}?access_token=${TOKEN}" >out_${FILENAME}.log

    echo
    nfile=$(( $nfile + 1 ))
done

echo Total number of files: $nfile
#--------------------
cd ${here}
