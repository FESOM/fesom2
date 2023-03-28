#!/usr/bin/env bash

# - - -
# # synopsis
# determine which environment directory is to be used for the current host
# # usage
# *source* to silently source the environment for this host system
# *execute* to print the environment directory for this host system 
# - - -


# see if we are being sourced or executed
# as we use bash to execute (see shebang), BASH_SOURCE is set when executing
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
   BEING_EXECUTED=true
else
   BEING_EXECUTED=false
fi

# if an arg is given and doesn't start with - use it as hostname, arguments with - are passed on to cmake
if [[ ! -z "$1" ]] && [[ ! "$1" = ^- ]]; then
   LOGINHOST=$1 # arg exists and doesn't start with -
   shift # pop the argument as we already stored it
else
   # no argument given
   LOGINHOST="$(hostname -f)"
fi

if [[ $LOGINHOST =~ ^m[A-Za-z0-9]+\.hpc\.dkrz\.de$ ]]; then
   STRATEGY="mistral.dkrz.de"
elif [[ $LOGINHOST =~ ^levante ]] || [[ $LOGINHOST =~ ^l[:alnum:]+\.lvt\.dkrz\.de$ ]]; then 
   STRATEGY="levante.dkrz.de"
   # following regex only matches if input is 2 word like levante.nvhpc, this enables using different shells for a machine directly
   compid_regex="^([[:alnum:]]+)\.([[:alnum:]]+)$"
   if [[ $LOGINHOST =~ $compid_regex ]]; then
     COMPILERID="${BASH_REMATCH[2]}"
   fi
elif [[ $LOGINHOST =~ ^ollie[0-9]$ ]] || [[ $LOGINHOST =~ ^prod-[0-9]{4}$ ]]; then
   STRATEGY="ollie"
elif [[ $LOGINHOST =~ ^albedo[0-9]$ ]] || [[ $LOGINHOST =~ ^prod-[0-9]{4}$ ]]; then
   STRATEGY="albedo"
elif [[ $LOGINHOST =~ ^h[A-Za-z0-9]+\.hsn\.hlrn\.de$ ]]; then
   STRATEGY="hlogin.hlrn.de"
elif [[ $LOGINHOST =~ ^b[A-Za-z0-9]+\.hsn\.hlrn\.de$ ]]; then
   STRATEGY="blogin.hlrn.de"
elif [[ $LOGINHOST =~ \.hww\.de$ ]] || [[ $LOGINHOST =~ ^nid[0-9]{5}$ ]]; then
   STRATEGY="hazelhen.hww.de"
elif [[  $LOGINHOST =~ \.jureca$ ]]; then
   STRATEGY="jureca"
elif [[  $LOGINHOST = ubuntu ]]; then
   STRATEGY="ubuntu"
elif [[  $LOGINHOST = bsc ]]; then
   STRATEGY="bsc"
elif [[  $LOGINHOST =~ ^juwels[0-9][0-9].ib.juwels.fzj.de$ ]]; then
   STRATEGY="juwels"
elif [[  $LOGINHOST =~ ^jwlogin[0-9][0-9].juwels$ ]]; then
   STRATEGY="juwels"
elif [[ $LOGINHOST =~ ^cc[a-b]+-login[0-9]+\.ecmwf\.int$ ]]; then
   STRATEGY="ecaccess.ecmwf.int"
elif [[ $LOGINHOST =~ ^stco-esl[0-9]+$ ]]; then
   STRATEGY="aleph"
elif [[ $LOGINHOST =~ ^[A-Za-z0-9]+\.ecmwf\.int$ ]]; then
STRATEGY="wsecmwf"
elif [[ $LOGINHOST =~ \.bullx$ ]]; then
STRATEGY="atosecmwf"
else
   echo "can not determine environment for host: "$LOGINHOST
   [ $BEING_EXECUTED = true ] && exit 1
   return # if we are being sourced, return from this script here
fi

if [ -n "$BASH_VERSION" ]; then
   # assume bash
   SOURCE="${BASH_SOURCE[0]}"
elif [ -n "$ZSH_VERSION" ]; then
   # assume zsh
   SOURCE=${(%):-%N}
fi

DIR="$( cd "$( dirname "${SOURCE}" )" && pwd )"

if [ $BEING_EXECUTED = true ]; then
   # file is being executed, why is this here?
   echo $DIR/env/$STRATEGY
else
   # file is being sourced
   export FESOM_PLATFORM_STRATEGY=$STRATEGY
   SHELLFILE="${DIR}/env/${STRATEGY}/shell"
   if [[ -n ${COMPILERID} ]]; then
      SHELLFILE="${SHELLFILE}.${COMPILERID}"
   fi
   if [[ ! -e ${SHELLFILE} ]]; then 
       echo "Shell file for ${LOGINHOST} doesnt exist: "$SHELLFILE
       exit 1
   fi
   source $SHELLFILE
fi
