#!/bin/bash
rm ./coawst*.out
rm ./COAWST_TIDE_OUTPUT/*
./predict_tide < ./COAWST_setup/setup_el_n.inp
./predict_tide < ./COAWST_setup/setup_u_n.inp
./predict_tide < ./COAWST_setup/setup_v_n.inp
./predict_tide < ./COAWST_setup/setup_el_e.inp
./predict_tide < ./COAWST_setup/setup_u_e.inp
./predict_tide < ./COAWST_setup/setup_v_e.inp
./predict_tide < ./COAWST_setup/setup_el_s.inp
./predict_tide < ./COAWST_setup/setup_u_s.inp
./predict_tide < ./COAWST_setup/setup_v_s.inp
#./predict_tide < ./COAWST_setup/setup_el_w.inp
#./predict_tide < ./COAWST_setup/setup_u_w.inp
#./predict_tide < ./COAWST_setup/setup_v_w.inp
./predict_tide < ./COAWST_setup/setup_el_ini.inp
./predict_tide < ./COAWST_setup/setup_u_ini.inp
./predict_tide < ./COAWST_setup/setup_v_ini.inp
mv ./coawst*.out ./COAWST_TIDE_OUTPUT
