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


case $choice1 in 1) 

# MENU 9 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Suspension data"

dialog --clear \
       --backtitle "$BACKTITLE" \
       --title "$TITLE" \
       --yesno "External magnetic field applied?" 7 60
       
       choice7=$?
clear
esac

case $choice1 in 3)  

# MENU 9 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="Suspension data"

dialog --clear \
       --backtitle "$BACKTITLE" \
       --title "$TITLE" \
       --yesno "External magnetic field applied?" 7 60
       
       choice7=$?
clear
esac

case $choice7 in 0)  
# MENU 10 

HEIGHT=15
WIDTH=60
CHOICE_HEIGHT=15
BACKTITLE="SIMMSUS - Simulation of Magnetic Suspensions"
TITLE="External field data"
MENU="Choose one of the following options for the applied field:"

OPTIONS=(1 "Steady-state field"
         2 "Oscillatory field"
         3 "Rotative field"
         4 "Oscillatory field with double frequency"
         5 "Nonlinear Duffing solution applied field"
         6 "Field-ramp (Dynamical increase)")

choice8=$(dialog --clear \
                --backtitle "$BACKTITLE" \
                --title "$TITLE" \
                --menu "$MENU" \
                $HEIGHT $WIDTH $CHOICE_HEIGHT \
                "${OPTIONS[@]}" \
                2>&1 >/dev/tty)
clear
esac

case $choice8 in
1) # MENU 11 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

alpha=""

exec 3>&1

form3=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Field data" \
--form "\nField data" 25 60 20 \
"Nondimensional field (alpha):" 1 1 "$alpha" 1 40 4 10 \
2>&1 1>&3)

exec 3>&-

alpha=$(echo "$form3" | sed -n 1p)
clear

esac

case $choice8 in
2) # MENU 12 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

alpha=""
omega_field=""

exec 3>&1

form4=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Field data" \
--form "\nField data" 40 100 25 \
"Nondimensional field (alpha):" 1 1 "$alpha" 1 40 4 25 \
"Nondimensional field frequency:" 2 1 "$omega_field" 2 40 4 25 \
2>&1 1>&3)

exec 3>&-

alpha=$(echo "$form4" | sed -n 1p)
omega_field=$(echo "$form4" | sed -n 2p)
clear

esac

case $choice8 in
3) # MENU 12 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

alpha=""
omega_field=""

exec 3>&1

form4=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Field data" \
--form "\nField data" 40 100 25 \
"Nondimensional field (alpha):" 1 1 "$alpha" 1 40 4 25 \
"Nondimensional field frequency:" 2 1 "$omega_field" 2 40 4 25 \
2>&1 1>&3)

exec 3>&-

alpha=$(echo "$form4" | sed -n 1p)
omega_field=$(echo "$form4" | sed -n 2p)
clear

esac

case $choice8 in
4) # MENU 13 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

alpha=""
omega_field=""
omega_field2=""

exec 3>&1

form5=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Field data" \
--form "\nField data" 40 100 25 \
"Nondimensional field (alpha):" 1 1 "$alpha" 1 40 4 25 \
"Nondimensional field frequency 1:" 2 1 "$omega_field" 2 40 4 25 \
"Nondimensional field frequency 2:" 3 1 "$omega_field2" 3 40 4 25 \
2>&1 1>&3)

exec 3>&-

alpha=$(echo "$form5" | sed -n 1p)
omega_field=$(echo "$form5" | sed -n 2p)
omega_field2=$(echo "$form5" | sed -n 3p)
clear

esac

case $choice8 in
5) # MENU 14 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

alpha=""
omega_field=""
c1_duffing=""
c2_duffing=""
c3_duffing=""
c4_duffing=""

exec 3>&1

form6=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Field data" \
--form "\nField data" 40 100 25 \
"Nondimensional field (alpha):" 1 1 "$alpha" 1 40 4 25 \
"Nondimensional field frequency:" 2 1 "$omega_field" 2 40 4 25 \
"C1 Duffing coefficient:" 3 1 "$c1_duffing" 3 40 4 25 \
"C2 Duffing coefficient:" 4 1 "$c2_duffing" 4 40 4 25 \
"C3 Duffing coefficient:" 5 1 "$c3_duffing" 5 40 4 25 \
"C4 Duffing coefficient:" 6 1 "$c4_duffing" 6 40 4 25 \
2>&1 1>&3)

exec 3>&-

alpha=$(echo "$form6" | sed -n 1p)
omega_field=$(echo "$form6" | sed -n 2p)
c1_duffing=$(echo "$form6" | sed -n 3p)
c2_duffing=$(echo "$form6" | sed -n 4p)
c3_duffing=$(echo "$form6" | sed -n 5p)
c4_duffing=$(echo "$form6" | sed -n 6p)
clear

esac

case $choice8 in
6) # MENU 15 
   
HEIGHT=15
WIDTH=40
CHOICE_HEIGHT=15

alpha=""
omega_field=""
max_omega=""
numberintdyn=""

exec 3>&1

form7=$(dialog --backtitle "SIMMSUS - Simulation of Magnetic Suspensions" --title "Field data" \
--form "\nField data" 40 100 25 \
"Nondimensional field (alpha):" 1 1 "$alpha" 1 40 4 25 \
"Nondimensional minimum frequency:" 2 1 "$omega_field" 2 40 4 25 \
"Nondimensional maximum frequency:" 3 1 "$max_omega" 3 40 4 25 \
"Number of intervals:" 4 1 "$numberintdyn" 4 40 4 25 \
2>&1 1>&3)

exec 3>&-

alpha=$(echo "$form7" | sed -n 1p)
omega_field=$(echo "$form7" | sed -n 2p)
max_omega=$(echo "$form7" | sed -n 3p)
numberintdyn=$(echo "$form7" | sed -n 4p)
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
    case $choice7 in
        0) echo "Choice 7 - option YES";; # some action on 1
        1) echo "Choice 7 - option NO";; # some action on 2
    esac
    case $choice8 in
        1) echo "Choice 8 - option 1 - alpha =" $alpha;; # some action on 1
        2) echo "Choice 8 - option 2 - alpha, omega =" $alpha, $omega_field;; # some action on 2
        3) echo "Choice 8 - option 3 - alpha, omega =" $alpha, $omega_field;; # some action on 2
        4) echo "Choice 8 - option 4 - alpha, omega1, omega2 =" $alpha, $omega_field, $omega_field2;; # some action on 2
        5) echo "Choice 8 - option 5 - alpha, omega, C1,C2,C3,C4 =" $alpha, $omega_field, $c1_duffing, $c2_duffing, $c3_duffing, $c4_duffing ;; # some action on 2
    	6) echo "Choice 8 - option 6 - alpha, omega_min, omega_max, number_intervals  =" $alpha, $omega_field, $max_omega, $numberintdyn ;; # some action on 2     
    esac
