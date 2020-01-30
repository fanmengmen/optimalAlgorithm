#include <iostream>
#include "inc/creatPlan.hpp"

using json = nlohmann::json;
using namespace std;



int PLAN::creatHead(const double hoverSpeed, const double cruiseSpeed, const double firmwareType)
{
    try
    {
        plan["fileType"] = "Plan";
        plan["geoFence"] = {{"circles",json::array()}, {"polygons",json::array()}, {"version",2}};
        plan["groundStation"] = "QGroundControl";
        plan["mission"] = {{"cruiseSpeed", cruiseSpeed}, {"firmwareType", firmwareType}, {"hoverSpeed",hoverSpeed}, {"items",{}}, {"plannedHomePosition",{}},{"vehicleType",2},{"version",2}};
        plan["rallyPoints"] = {{"points", json::array()},{"version",2}};
        plan["version"] = 1;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return -1;
    }
    
    
    return 0;
}


int PLAN::setHomePosition(const double home_lat, const double home_lon, const double home_alt)
{
    plan["mission"]["plannedHomePosition"] = {home_lat, home_lon, home_alt};
}

bool PLAN::addWayPoint(int num, double lat, double lon, double alt)
{
    json tmp;
    tmp["AMSLAltAboveTerrain"] = {};
    tmp["Altitude"] = alt;
    tmp["AltitudeMode"] = 1;
    tmp["autoContinue"] = true;
    tmp["command"] = 16;
    tmp["doJumpId"] = num;
    tmp["frame"] = 3;
    tmp["params"] = {0, 0, 0, {}, lat, lon, alt};
    tmp["type"] = "SimpleItem";
    // cout << setw(4) << tmp << endl;
    plan["mission"]["items"].push_back(tmp);
}


void  PLAN::printJson(){
    cout << setw(4) << plan << endl;
}


int PLAN::outJson(const string fileName)
{
    
    ofstream out;
    out.open(fileName);
    out << setw(4) << plan << endl;
    out.close();
   
    return 0;
    
}