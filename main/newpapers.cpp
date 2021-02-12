#include "Cia402device.h"
#include "CiA301CommPort.h"
#include "SocketCanPort.h"
#include <iostream>
#include <stdio.h>

#include <fstream>
#include <chrono>
#include <ctime>
#include <ios>

#include <boost/asio.hpp> // include boost
#include <boost/asio/serial_port.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <string.h>
#include <math.h>
#include <sstream>
#include <boost/algorithm/hex.hpp>
#include "imu3dmgx510.h"
#include "OnlineSystemIdentification.h"

#include "math.h"

#include "fcontrol.h"
#include "IPlot.h"
#include "Kinematics.h"


#include <tuple>

using namespace boost::asio;
using namespace boost::algorithm;
using namespace std::string_literals;
using namespace stateestimation;
using std::cin;
using std::cout;

#ifdef _WIN32
// windows uses com ports, this depends on what com port your cable is plugged in to.
const char *PORT = "COM7";
#else
//default port usb
const char *PORT = "/dev/ttyUSB0";
#endif

void setup(){

    SocketCanPort pm31("can1");
    CiA402SetupData sd31(2048,24,0.001, 0.144, 20);
    CiA402Device m1 (1, &pm31, &sd31);
    m1.Reset();
    m1.SwitchOn();
    m1.SetupPositionMode(10,10);

    SocketCanPort pm2("can1");
    CiA402SetupData sd32(2048,24,0.001, 0.144, 20);
    CiA402Device m2 (2, &pm2, &sd32);
    m2.Reset();
    m2.SwitchOn();
    m2.SetupPositionMode(10,10);

    SocketCanPort pm3("can1");
    CiA402SetupData sd33(2048,24,0.001, 0.144, 20);
    CiA402Device m3 (3, &pm3, &sd33);
    m3.Reset();
    m3.SwitchOn();
    m3.SetupPositionMode(10,10);

    IMU3DMGX510 misensor("/dev/ttyUSB0");

    misensor.set_IDLEmode();
    misensor.set_freq(10);
    misensor.set_devicetogetgyroacc();
    misensor.set_streamon();
    cout << "Calibrating IMU..." << endl;
    misensor.calibrate();
    cout << "Calibration done" << endl;
    double *EulerAngles;

    misensor.set_streamon();
    int cnt=0;

    do{

        EulerAngles = misensor.EulerAngles();

        cout << "Roll: " << EulerAngles[0] << " Pitch: " << EulerAngles[1] << endl;
        cnt=cnt+1111;

    }while(cnt>0);



    double pos;
    double vel;
    cout << m1.SetPosition(0) << endl;
    cout << m2.SetPosition(0) << endl;
    cout << m3.SetPosition(0) << endl;

    // position  [rads]
    cout << m1.GetPosition() << endl;
    cout << m2.GetPosition() << endl;
    cout << m3.GetPosition() << endl;
}

void capturedata(){
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    time_t tt = std::chrono::system_clock::to_time_t(now);
    tm local_tm = *localtime(&tt);
    cout<<local_tm.tm_mon<<" "<<local_tm.tm_mday<<" "<<local_tm.tm_hour<<endl;
    string address="/home/humasoft/code/papers/graficas/newpaper/Dataset"+to_string(local_tm.tm_mon)+"_"+to_string(local_tm.tm_mday)+"_"+to_string(local_tm.tm_hour)+"_"+to_string(local_tm.tm_min)+".csv";
    ofstream data(address,std::ofstream::out);

    SocketCanPort pm31("can1");
    CiA402SetupData sd31(2048,24,0.001, 0.144, 20);
    CiA402Device m1 (1, &pm31, &sd31);
    m1.Reset();
    m1.SwitchOn();
    m1.SetupPositionMode(10,10);

    SocketCanPort pm2("can1");
    CiA402SetupData sd32(2048,24,0.001, 0.144, 20);
    CiA402Device m2 (2, &pm2, &sd32);
    m2.Reset();
    m2.SwitchOn();
    m2.SetupPositionMode(10,10);

    SocketCanPort pm3("can1");
    CiA402SetupData sd33(2048,24,0.001, 0.144, 20);
    CiA402Device m3 (3, &pm3, &sd33);
    m3.Reset();
    m3.SwitchOn();
    m3.SetupPositionMode(10,10);

    IMU3DMGX510 misensor("/dev/ttyUSB0");

    misensor.set_IDLEmode();
    misensor.set_freq(10);
    misensor.set_devicetogetgyroacc();
    misensor.set_streamon();
    cout << "Calibrating IMU..." << endl;
    misensor.calibrate();
    cout << "Calibration done" << endl;
    double *EulerAngles;

    misensor.set_streamon();
    double dts=0.02;
    double f=0;

    SamplingTime Ts(dts);

    double psr = 0.0, isignal = 0.0;

    for(double t=dts;t<100;t=t+dts)
    {
        f=f+0.002;
        isignal = (1.5+3*sin(f*t));
        m3.SetPosition(isignal);
        //m2.SetPosition(0);
        cout << "t: "<< t << ", pos: " << isignal << endl;
        Ts.WaitSamplingTime();


        cout<<"Read position: "<<m3.GetPosition()<<", vel: "<<m3.GetVelocity()
            <<" and those amps:"<<m3.GetAmps()<<endl;

        EulerAngles = misensor.EulerAngles();

        cout << "ROLL: " << EulerAngles[0] << " ; PITCH: "  << EulerAngles[1] << endl;

        data << t << ", "  << isignal << ", "<< m3.GetPosition() <<", "<< m3.GetVelocity()
             <<", "<< m3.GetAmps() <<", "<<  EulerAngles[0] << ", " << EulerAngles[2] << endl;
    }


}
int main(){
    capturedata();
    setup();
}
