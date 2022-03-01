#!/bin/bash
## Create a gif of individually saved frames and label the frame number

## Assumes files named ${FNAME}${FRAME}.png
FNAME=my_vmd_image
START_FRAME=0
END_FRAME=12

## ${OUT_GIF}.${EXT}
## Pre-extension name for the output gif
OUT_GIF=my_output_file
## Lowercase extension
EXT=tga

#---------------------------- Behind the Curtain ----------------------------#
## Set up variables for loop
CURR_FRAME=${START_FRAME}
W_END_FRAME=$[${END_FRAME}+1]

while [ ${CURR_FRAME} -lt ${W_END_FRAME} ]; do
  ## Add the frame number to the top right corner, first in black for outline,
  ##  then in white. Resave as a new image.
  ## VMD TGA files get saved weirdly and need to be flipped
  if [[ "$EXT" == "tga" ]]
  then
    convert ${FNAME}${CURR_FRAME}.${EXT} -gravity northeast \
      -font Arial -stroke '#000C' -strokewidth 5 -pointsize 48 \
         -annotate 0 ${CURR_FRAME} \
      -font Arial -stroke  none   -fill white    -pointsize 48 \
         -annotate 0 ${CURR_FRAME} -flip \
      ${FNAME}.${CURR_FRAME}.${EXT}
  else
  ## For PNGs
    convert ${FNAME}${CURR_FRAME}.${EXT} -gravity northeast \
      -font Arial -stroke '#000C' -strokewidth 5 -pointsize 48 \
         -annotate 0 ${CURR_FRAME} \
      -font Arial -stroke  none   -fill white    -pointsize 48 \
         -annotate 0 ${CURR_FRAME} \
      ${FNAME}.${CURR_FRAME}.${EXT}
  fi
  ## Increase counter
  CURR_FRAME=$[${CURR_FRAME}+1]
done

## Get File List for gif (in order)
## Creates the list of files named ${FNAME}.${FRAME}.png
FILES=$(seq -f "${FNAME}.%g.${EXT}" -s " " ${START_FRAME} ${END_FRAME})

## Convert into a continuously looping gif (0), pausing for 30 ticks
convert -delay 30 -layers optimize -dispose previous -loop 0 \
  ${FILES} \
  ${OUT_GIF}.gif
