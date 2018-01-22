# nilmtk-lvns
The compatibility between Non-Intrusive Load Monitoring Toolkit and Low Voltage Network Simulator: A fast way to improve development of energy disaggregation algorithms

One obstacle in developing Non-Intrusive Load Monitoring (NILM) techniques is the difficulty of getting good home data set. A good data set consists of simultaneous measurements of the main branch and the appliances present in the residence so that the estimation obtained by NILM techniques from the main branch can be compared with the values measured for each appliance. Because it is necessary to install a meter on the main branch of the residence and a meter on each appliance, this is an expensive, laborious task and is often subject to noise and inaccuracies. The Non-Intrusive Load Monitoring Toolkit has made good progress by proposing a common data set format into which existing publicly available data sets can be converted, thereby enabling comparison of energy disaggregation algorithms in a reproducible manner. In a complementary way, the Low Voltage Network Simulator generates a simulation of one day of several residences, allowing an evaluation of several scenarios and avoiding problems like noise and lack of synchrony of the measurements. This paper proposes the compatibility of the both tools in order to facilitate the development and testing of new NILM techniques. The compatibility allows an agile evaluation of several scenarios and the generation of ground truth for the intermediate stages of the NILM, thus enabling a precise evaluation that avoids the propagation of errors between the stages

# Install
* Install [Non-Intrusive Load Monitoring Toolkit](https://github.com/nilmtk/nilmtk/blob/master/docs/manual/user_guide/install_user.md) 
* Obtain the [Low Voltage Network Simulator](http://www.dsee.fee.unicamp.br/~torquato/) and extract.
* Download the files `main_home_simulator.py` and `definitions.py`. Replace them in the `Low_Voltage_Network_Simulator_v1_Python` folder.
* After place the folder `lvns` in folder of converters of NILMTK, in general located in `Anaconda\envs\nilmtk-env\Lib\site-packages\nilmtk\dataset_converters`

After finish you can run `Convert_LVNS_3houses_base_to_NILMTK.ipynb` as example
