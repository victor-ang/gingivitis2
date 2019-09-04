#!/bin/bash
QT_INSTALL_DIR=/home/cussat/Qt/5.12.2/gcc_64
export Qt5_DIR=${QT_INSTALL_DIR}/lib/cmake
export Qt5Core_DIR=${QT_INSTALL_DIR}/lib/cmake/Qt5Core
export Qt5Widgets_DIR=${QT_INSTALL_DIR}/lib/cmake/Qt5Widgets
export Qt5Test_DIR=${QT_INSTALL_DIR}/lib/cmake/Qt5Test
export Qt5Gui_DIR=${QT_INSTALL_DIR}/lib/cmake/Qt5Gui
export Qt5Concurrent_DIR=${QT_INSTALL_DIR}/lib/cmake/Qt5Concurrent
#export QT_DEBUG_PLUGINS=1

mkdir -p build; cd build;

main(){
    echo -e "\e[34m ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ "
    echo -e "\e[34m ⟢  \e[34m        Builder   \e[34m ⟢  "
    echo -e "\e[34m ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ "
    echo ""
    echo -e "\e[0m 1\e[31m ▹\e[0m Compile"
    echo -e "\e[0m 2\e[31m ▹\e[0m Clean Build Compile"
    echo -e "\e[0m 3\e[31m ▹\e[0m Run"
    echo -e "\e[0m 4\e[31m ▹\e[0m Commit Changes"
    echo -e "\e[0m 5\e[31m ▹\e[0m Push last version"
    echo -e "\e[0m 6\e[31m ▹\e[0m Exit"
    
    echo -e "\e[96m ?"
    while : ; do
        read input
        if [ \( "$input" -lt 7 \) -a \( "$input" -ge 1 \) ]; then
            break
        else
            echo -e "\e[31mout of range\e[0m"
        fi
    done
    echo -e "\e[0m "

    case $input in
        1 ) clear; clear; make -j clean; make -j ;  ;;
		2 ) clear; clear; rm -rf *; cmake ..; make -j ;  ;;
        3 ) ./viewer;;
        4 ) git add ../*; git commit;;
        5 ) git push ;;
        6 ) exit;;
    esac
}

launcher(){
    clear
    echo -e "\e[34m ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ "
    echo -e "\e[34m ⟢  \e[34m     Launcher     \e[34m ⟢  "
    echo -e "\e[34m ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ ⟢ "
    list="$(find . -maxdepth 1 -type f -perm +111)"
    declare -i i=1
    declare -a cmd=()
    IFS=$'\n' sorted=($(sort <<<"${list[*]}"))
    unset IFS
    for file in ${sorted[@]}; do 
        echo -e "\e[0m $i\e[31m ▹\e[0m ${file//.\/}"
        cmd+=($file)
        i+=1
        if [[ $file = *"viewer"* ]]; then
            echo -e "\e[0m $i\e[31m ▹\e[0m ${file//.\/} \e[31m[WITH GPU]"
            cmd+=("optirun $file")
            i+=1
        fi
    done
    echo -e "\e[96m ?"
    while : ; do
        read input
        if [ \( "$input" -lt "$i" \) -a \( "$input" -ge 1 \) ]; then
            break
        else
            echo -e "\e[31mout of range\e[0m"
        fi
    done

    clear
    echo -e "\e[34m ⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢  "
    echo -e "\e[34m   \e[34m ${cmd[$input-1]}"
    echo -e "\e[34m ⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢⟢  "
    echo -e "\e[0m "
    exec ${cmd[$input-1]}
}

clear
while : ; do
    main
done
