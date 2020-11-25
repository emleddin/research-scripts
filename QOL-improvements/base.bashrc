
## Check the window size after each command and update the values of 
## LINES and COLUMNS if necessary.
shopt -s checkwinsize

## Bashrc aliases
alias prof='vi ~/.bashrc'
alias sorcery='source ~/.bashrc && echo "Your ~/.bashrc has been sourced, sire."'

alias ana="python3.6"
alias n="nvidia-smi"
alias dirt="tree -d"
alias lj="ls -lthr"
alias lf="ls -l1vh"
alias tails="tail -n 15"
alias mymail="cat /var/mail/username"

## PBS Queue Scheduler
alias wd="cd $PBS_O_WORKDIR"
alias qme="qstat -u $USER -n"
alias jobs="qstat -u $USER"
alias sehen="watch -n 1 qstat -u $USER -n"

## SLURM Queue Scheduler
alias wd="cd $SLURM_SUBMIT_DIR"
alias showq='squeue -u $USER'
alias qme='squeue -u $USER'
alias sehen='watch -n 5 squeue -u $USER'

## Storage nonsense
alias storespace="df -h"
alias storeuse="du -sh *"

## Make images
alias eps2png='for file in *.eps; do convert $file -rotate 90 ${file%.*}.png; done'

## Typo Aliases
alias sheen="watch -n 1 qstat -u $USER"
alias pdw="pwd"
alias tial="tail"

## Appearance Aliases
## Change the parts in <> if you only need to change name sometimes
if [ "$HOSTNAME" = <specific.login.hostname> ]; then
    HN="<nicename>"
else
    HN=$(hostname -s)
fi
## If you have a specific name you always want to use, use this instead
HN="<use_this_name>"
## Set up the prompt appearance
PS1='\[\e[0;32m\]\u@${HN}\[\e[0;37m\]:\[\e[0;36m\]\w\[\e[0;37m\]$ '
## No color, but full path
#PS1='\u@${HN}:\w$ '
## No color and last folder only
#PS1='\u@\h:\W$ '

## Depending on installs and globals
alias vi='vim'

