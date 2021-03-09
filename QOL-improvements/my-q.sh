get_running () {
  # Set the HN (change a non-descript name to something meaningful)
  # Or use the descriptive one in short form
  if [ "$HOSTNAME" = long-nondescript-hostname ]; then
      HN="computerA"
  else
      HN=$(hostname -s)
  fi

  printf "           ${USER}'s jobs on ${HN}\n\n"

  # Print the Header
  printf "%-7s %-8s %-9s %-18s %-8s %-2s\n" \
  "JOBID" "USER" "QUEUE" "Job Name" "Node" "Status"
  echo "================================================================"

  # Get and print the running
  qstat -u $USER -n | grep -A1 -w R | egrep -v "^ " | \
  awk '{sub(/\.cruntch3.chem.u.*$/,"",$1); \
  printf "%-7s %-8s %-9s %-18s %-8s %-2s\n", \
  $1,$2,$3,$4,$7,$10}'

  # Get and print the queued
  qstat -u $USER -n | grep -A1 -w Q | egrep -v "^ " | \
  awk '{sub(/\.cruntch3.chem.u.*$/,"",$1); \
  printf "%-7s %-8s %-9s %-18s %-8s %-2s\n", \
  $1,$2,$3,$4,$7,$10}'
}

get_running
