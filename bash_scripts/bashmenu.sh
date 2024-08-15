#!/bin/bash

# MENU 1 

HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Particle data"
MENU="Choose one of the following options:"

OPTIONS=(1 "Magnetic particles"
         2 "Non-magnetic particles"
         3 "A mixture of both")

choice1=$(dialog --clear \
                --backtitle "$BACKTITLE" \
                --title "$TITLE" \
                --menu "$MENU" \
                $HEIGHT $WIDTH $CHOICE_HEIGHT \
                "${OPTIONS[@]}" \
                2>&1 >/dev/tty)
clear

# MENU 2 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Particle data"
MENU="Choose one of the following options:"

OPTIONS=(1 "Brownian particles"
         2 "Non-Brownian particles")

choice2=$(dialog --clear \
                --backtitle "$BACKTITLE" \
                --title "$TITLE" \
                --menu "$MENU" \
                $HEIGHT $WIDTH $CHOICE_HEIGHT \
                "${OPTIONS[@]}" \
                2>&1 >/dev/tty)
clear

# MENU 3 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Particle data"
MENU="Choose one of the following options:"

OPTIONS=(1 "Massive particles"
         2 "Non-massive particles")

choice3=$(dialog --clear \
                --backtitle "$BACKTITLE" \
                --title "$TITLE" \
                --menu "$MENU" \
                $HEIGHT $WIDTH $CHOICE_HEIGHT \
                "${OPTIONS[@]}" \
                2>&1 >/dev/tty)
clear

# MENU 4 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Particle data"
MENU="Choose one of the following options:"

OPTIONS=(1 "Neutrally buoyant particles"
         2 "Sedimenting particles")

choice4=$(dialog --clear \
                --backtitle "$BACKTITLE" \
                --title "$TITLE" \
                --menu "$MENU" \
                $HEIGHT $WIDTH $CHOICE_HEIGHT \
                "${OPTIONS[@]}" \
                2>&1 >/dev/tty)
clear

# MENU 5 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Suspension data"

dialog --clear \
       --backtitle "$BACKTITLE" \
       --title "$TITLE" \
       --yesno "Suspension subjected to shear-rate?" 7 60
       
       choice5=$?
clear

    case $choice5 in
0) # MENU 6 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Suspension data"

dialog --clear \
       --backtitle "$BACKTITLE" \
       --title "$TITLE" \
       --yesno "Oscillatory shear-imposed?" 7 60
       
       choice6=$?
clear
esac

case $choice6 in
0) # MENU 7 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

shearamplitude=""
shearfrequency=""

exec 3>&1

form1=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Suspension data" \
--form "\nShear rate data" 25 60 16 \
"Shear rate amplitude:" 1 1 "$shearamplitude" 1 25 25 30 \
"Shear rate frequency:" 2 1 "$shearfrequency" 2 25 25 30 \
2>&1 1>&3)

exec 3>&-

shearamplitude=$(echo "$form1" | sed -n 1p)
shearfrequency=$(echo "$form1" | sed -n 2p)
clear

esac

case $choice6 in
1) # MENU 8 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

shearamplitude=""

exec 3>&1

form2=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Suspension data" \
--form "\nShear rate data" 25 60 16 \
"Shear rate amplitude:" 1 1 "$shearamplitude" 1 25 25 30 \
2>&1 1>&3)

exec 3>&-

shearamplitude=$(echo "$form2" | sed -n 1p)
clear

esac



   
    case $choice1 in
        1) echo "Choice 1 - option 1";; # some action on 1
        2) echo "Choice 1 - option 2";; # some action on 2
        3) echo "Choice 1 - option 3";; # some action on other
    esac
    
    case $choice2 in
        1) echo "Choice 2 - option 1";; # some action on 1
        2) echo "Choice 2 - option 2";; # some action on 2
    esac
    
    case $choice3 in
        1) echo "Choice 3 - option 1";; # some action on 1
        2) echo "Choice 3 - option 2";; # some action on 2
    esac
    
    case $choice4 in
        1) echo "Choice 4 - option 1";; # some action on 1
        2) echo "Choice 4 - option 2";; # some action on 2
    esac
    case $choice5 in
        0) echo "Choice 5 - option YES";; # some action on 1
        1) echo "Choice 5 - option NO";; # some action on 2
    esac
    case $choice6 in
        0) echo "Shear amplitude and frequency:" $shearamplitude, $shearfrequency;; # some action on 1
        1) echo "Shear amplitude:" $shearamplitude;; # some action on 1
    esac
