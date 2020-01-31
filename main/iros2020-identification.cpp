

#include <iostream>
#include <math.h>

#include "TestPort.h"
#include "CanBusPort.h"
#include "Cia402device.h"
#include "CiA301CommPort.h"
#include "SocketCanPort.h"
#include "fcontrol.h"
#include "SensorIntegration.h"



using namespace std;

int main(int argc, char *argv[])
{

    ///prepare ports
    /// Open a port address with a PortBase Object
    //TestPort port;
    //CanBusPort port;
    SocketCanPort p1("can1");


    ///Then instantiate a 301 communications object for that address
    //CiA301CommPort coms(port.getPortFileDescriptor());

    ///Create a joint and give a canopen id, and a 301port (by constructor)
    CiA402Device m1(32,&p1);

    ///--sensor tilt--
    SerialArduino tilt;
    double incSensor,oriSensor;


    if (tilt.readSensor(incSensor,oriSensor) <0){}
    while(true)
        cout << "Inc: " << incSensor << " ; Ori: "  << oriSensor << endl;


    m1.SetupPositionMode(200,200);
//    m1.Setup_Velocity_Mode(0,360);
//    m1.Setup_Torque_Mode();


    double dts=0.01;
    SamplingTime Ts(dts);



    for(double t=dts;t<10;t=t+dts)
    {
        m1.SetPosition(3+3*sin(5*t));
        cout << "t: "<< t << ", pos: " << +3*sin(5*t) << endl;
        Ts.WaitSamplingTime();

    }


    return 0;

}
