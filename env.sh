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

# if an arg is given, use it as hostname
if [ -z "$1" ]; then
   # no argument given
   LOGINHOST="$(hostname -f)"
else
   LOGINHOST=$1
fi

if [[ $LOGINHOST =~ ^m[A-Za-z0-9]+\.hpc\.dkrz\.de$ ]]; then
   STRATEGY="mistral.dkrz.de"
elif [[ $LOGINHOST =~ ^ollie[0-9]$ ]] || [[ $LOGINHOST =~ ^prod-[0-9]{4}$ ]]; then
   STRATEGY="ollie"
elif [[ $LOGINHOST =~ ^h[A-Za-z0-9]+\.hsn\.hlrn\.de$ ]]; then
   STRATEGY="hlogin.hlrn.de"
elif [[ $LOGINHOST =~ ^b[A-Za-z0-9]+\.usr\.hlrn\.de$ ]]; then
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
   # file is being executed
   echo $DIR/env/$STRATEGY
else
   # file is being sourced
   source $DIR/env/$STRATEGY/shell
fi
