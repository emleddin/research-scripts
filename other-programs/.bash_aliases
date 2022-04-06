# This file contains aliases to add to a ~/.bashrc file
# or to add a command to source.

#---- Internal directories
#alias rn='cd '

#----- ssh
## XSEDE
#alias xsede='ssh user@login.xsede.org'

#----- Keyboard Shortcuts
alias lj="ls -lthr"
alias lf="ls -l1vh"
alias tails="tail -n 15"
alias n="nvidia-smi"
alias dirt="tree -d"
alias fopen='nautilus .'
alias desk='cd ~/Desktop'
alias down='cd ~/Downloads'
alias storespace="df -h"
alias storeuse="du -sh *"

## Computer Management
alias u1='sudo apt update && sudo apt upgrade'
alias u2='sudo apt autoremove'
alias sshres='sudo service ssh restart'

#----- Typo Aliases
alias pdw="pwd"
alias tial="tail"
alias mdkir="mkdir"
alias mkidr="mkdir"
alias ego="eog"

#----- Bashrc aliases
alias prof='vi ~/.bashrc'
alias sorcery='source ~/.bashrc && echo "Your ~/.bashrc was sourced, sire."'

#----- Make images
alias eps2png='for file in *.eps; do convert $file -rotate 90 ${file%.*}.png; done'
#alias eps2png='for file in *.eps; do convert $file ${file%.*}.png; done'

## Fix PPM to POV
alias ppm2png='for file in *.ppm; do mv $file ${file%.*}.png; done'

## Fix POV to PNG
alias pov2png='for filename in *pov.png; do mv $filename ${filename//.pov.png/.png}; done'

#----- Functions
## SFTP
sftp_system () {
sftp user@system:${1}
}

##Show Amber Averages
amber_avg () {
   grep -A 12 "A V E R A G E S" $1
}

## Convert au to kcal/mol
au2kcal () {
    echo "${1} * 627.51" | bc -l
}

## Products - Reactants
mysub () {
    echo "${1} - ${2}" | bc -l
}

## Ghostscript (`comp_pdf input output`)
comp_pdf () {
  gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dNOPAUSE -dQUIET -dBATCH \
  -sOutputFile=${2} ${1}
}
