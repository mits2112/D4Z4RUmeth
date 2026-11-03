
check_mkdir()
{
if [ -d $1 ]
then
    echo "$1 exists."
else
    echo "$1 does not exits."
    mkdir -p $1
fi
}

check_error()
{
if [ $1 -ne 0 ]; then
    echo "FATAL ERROR: pipeline script"
    echo "ERROR CODE: $1"
    exit $1
fi
}

check_file_exists()
{
if [ -f $1 ]; then
  echo "$1 exists."
  return
fi
echo "$1 does not exists."
exit 1
}




# Function to check job dependencies and set the 'dep' variable
check_dependency() {
    local job_name="$1"  # Argument passed to the function (job name)
    local job_id=""
    local dependency_arg=""

    if [ -n "$job_name" ]; then
        # Get the job ID using the squeue command
        # awk 'NR>1 {print $1}' extracts the first column ($1, the JOBID) excluding the header row (NR>1)
        job_id=$(squeue --name="$job_name" -h | awk '{print $1}')
        
        if [ -n "$job_id" ]; then
            # If a job ID is obtained, set the dependency argument
            dependency_arg="--dependency=afterok:$job_id"
            echo "Dependency set: $dependency_arg (Job name: $job_name, ID: $job_id)"
        else
            # If the job is not found
            dependency_arg=""
            echo "Warning: The specified job name '$job_name' was not found in the queue."
        fi
    else
        # If the job name is empty
        dependency_arg=""
        echo "prejob is not specified. No dependency."
    fi

    # Set the value from within the function to the global 'dep' variable (to make it available outside the function)
    dep="$dependency_arg"

}




