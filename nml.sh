#!/usr/bin/env bash

# dependencies
SOURCE_FILES=(nml.c nml_util.c)
OBJECT_FILES=(nml.o nml_util.o)
HEADER_FILES=(nml.h nml_util.h)
MAKEFILE=Makefile

# commands
CMDS=($CC ar nm)

# colored print
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# exit when any command fails
set -e

function sanity_checks {
  echo -e "${YELLOW}Sanity Checks:${NC}"

  # Check for existence of required commands
  for cmd in ${CMDS[@]};
  do
    echo -ne "\tChecking if $cmd is in PATH:"
    if ! command -v ${cmd} &> /dev/null
    then
      echo -e " ${RED}NOK${NC}"
      echo -e " ${RED}ABORTING${NC}"
      exit
    fi
    echo -e " ${GREEN}OK${NC}"
  done

  # check for source files
  for file in ${SOURCE_FILES[@]};
  do
    echo -ne "\tChecking if src/${file} exists in folder:"
    if [ ! -f "src/${file}" ]; then
      echo -e " ${RED}NOK${NC}"
      echo -e " ${RED}ABORTING${NC}"
      exit
    fi
    echo -e " ${GREEN}OK${NC}"
  done

  # check for header files
  for header in ${HEADER_FILES[@]};
  do
    echo -ne "\tChecking if include/${header} exists in folder:"
    if [ ! -f "include/${header}" ]; then
      echo -e " ${RED}NOK${NC}"
      echo -e " ${RED}ABORTING${NC}"
      exit
    fi
    echo -e " ${GREEN}OK${NC}"
  done

  # check for Makefile
  test -e ${MAKEFILE}
}

function usage {
  echo -e "Usage:"
  echo -e "\t ${YELLOW}./nml.sh build${NC}"
  echo -e "\t\t Builds the lib in: ${YELLOW}${DIST_DIR}/${NC}."
  echo -e "\t ${YELLOW}./nml.sh tests${NC}"
  echo -e "\t\t Runs lib tests."
  echo -e "\t ${YELLOW}./nml.sh examples${NC}"
  echo -e "\t\t Builds ${YELLOW}${EXAMPLES}${NC}/ folder with the latest build."
  echo -e "\t ${YELLOW}./nml.sh clean${NC}"
  echo -e "\t\t Cleans the folder for *.o and *.a files. Deletes the ${YELLOW}${DIST_DIR}/${NC} folder."
}

### MAIN ###

echo -e " "
echo -e "${GREEN}N${NC}eat ${GREEN}M${NC}atrix ${GREEN}L${NC}ibrary (libnml)"
echo -e " "

for token in "${@}" ;
do
  case $token in
  "check")
    sanity_checks
    ;;
  "build")
    sanity_checks
    make objects lib dist
    ;;
  "examples")
    sanity_checks
    make examples
    ;;
  "test")
    sanity_checks
    make tests
    ;;
  "clean")
    make clean
    ;;
  "all")
    make all
    ;;
  *)
    echo -e "${RED}Unknown Option: '${1}'.${NC}"
    usage
    ;;
  esac
done
