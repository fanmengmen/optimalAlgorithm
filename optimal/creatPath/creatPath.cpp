#include<iostream>
#include<math.h>
#include<vector>
#include <fstream>
#include"inc/creatPlan.hpp"
#define pi 3.1415926
#define C_EARTH (double)6378137.0
#define DEG2RAD 0.0174532925
using namespace std;

int main()
{

    double speed = 3;
    double r = 15;
    // 中心点位置
    double center_lat = 30.7461401 * DEG2RAD;
    double center_lon = 103.9292323 * DEG2RAD;
    double center_alt = 530;
    
    // 起始点位置
    double origin_lat = 30.7457830 * DEG2RAD;
    double origin_lon = 103.9290763 * DEG2RAD;
    double origin_alt = 530;

    // 根据DJI的公式和数据
    // deltaNed.x      = deltaLat * C_EARTH;
    // deltaNed.y      = deltaLon * C_EARTH * cos(broadcastTarget->latitude);
    // deltaNed.z      = broadcastTarget->altitude - broadcastOrigin->altitude;
    


    // 根据pixhawk中的map_projection_project函数，
    /*1. 先将GPS的单位°转化为弧度（rad），因为c++中三角函数都是用的弧度
      2. 等距方位投影
         arg = ref->sin(lat) * sin(lat) + ref->cos(lat) * cos(lat) * cos(lon - ref->lon)
         if(arg > 1.0) arg = 1.0;
         else if(arg < -1)* arg = -1.0
         c = acos(arg)
         k = fabs(c) < DBL_EPSILON ? 1.0 : (c/sin(c))
         x = k*(cos(ref->lat)*sin(lat) - sin(ref->lat)*cos(lat)*cos(lon - ref->lon))*RADIUS_OF_EARTH;
         y = k*(cos(lat)*sin(lon - ref->lon)* RADIUS_OF_EARTH);
         */

    // double c = (center_lat - origin_lat) * C_EARTH;
    // double d = (center_lon - origin_lon) * C_EARTH * cos(center_lat);
    // cout << "c: " << c << "\t" << "d: " << d << endl;
    double c = 50;
    double d = 20;


    double sampleTime = 1;
    vector<double> data_x;
    vector<double> data_y;
    
    double T = 2*r/speed;
    //No.1 直线加圆弧轨迹构建
    for(double t = 0; t < T; t=t + sampleTime)
    {
      double x = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
      double y = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
      cout << x << "  " << y << endl;
      data_x.push_back((x + c)/C_EARTH + origin_lat);
      data_y.push_back((y + d)/C_EARTH/cos((x + c)/C_EARTH + origin_lat) + origin_lon);
    }

    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t + sampleTime)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      double x = sqrt(2)*r + r*cos(theta1);
      double y = r*sin(theta1);
      cout << x << "  " << y << endl;
      data_x.push_back((x + c)/C_EARTH + origin_lat);
      data_y.push_back((y + d)/C_EARTH/cos((x + c)/C_EARTH + origin_lat) + origin_lon);
    }

    T = 2*r/speed;
     for(double t = 0; t < T; t=t + sampleTime)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      double x = sqrt(2)/2*r - sqrt(2)/2*speed*t;
      double y = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
      cout << x << "  " << y << endl;
      data_x.push_back((x + c)/C_EARTH + origin_lat);
      data_y.push_back((y + d)/C_EARTH/cos((x + c)/C_EARTH + origin_lat) + origin_lon);
      
    }

    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t + sampleTime)
    {
      double theta2 = pi/4 + speed*t/r;
      double x = -sqrt(2)*r + r*cos(theta2);
      double y = r*sin(theta2);
      cout << x << "  " << y << endl;
      data_x.push_back((x + c)/C_EARTH + origin_lat);
      data_y.push_back((y + d)/C_EARTH/cos((x + c)/C_EARTH + origin_lat) + origin_lon);     
    }

    /*老版本的生成plan文件是通过字符串输出的方式，代码量太大*/
    // ofstream plan("eight.plan");  
    // plan.setf(ios::fixed, ios::floatfield);  // 设定为 fixed 模式，以小数点表示浮点数
    // plan.precision(14);  // 设置精度 2
    // if (plan.fail())  
    // {  
    //   cout<< "打开文件错误!" <<endl;  
    //   exit(0);  
    // }  
    // string str0 = "{\n\t\"fileType\": \"Plan\",\n\t\"geoFence\": {\n\t\t\"circles\": [\n\t\t],\n\t\t\"polygons\": [\n\t\t],\n\t\t\"version\": 2\n\t},\n\t\"groundStation\": \"QGroundControl\",\n\t\"mission\": {\n\t\t\"cruiseSpeed\": 15,\n\t\t\"firmwareType\": 12,\n\t\t\"hoverSpeed\": 3,\n\t\t\"items\": [\n";
    // string str1 = "{\n\t\"AMSLAltAboveTerrain\": null,\n\t\"Altitude\": 15,\n\t\"AltitudeMode\": 1,\n\t\"autoContinue\": true,\n\t\"command\": 16,\n\t\"doJumpId\": ";
    // string str2 = "\t\"frame\": 3,\n\t\"params\": [0, 0, 0, null, ";
    // string str3 = "{\n\t\"autoContinue\": true,\n\t\"command\": 178,\n\t\"doJupmId\": ";
    // string str4 = "\t\"frame\": 2,\n\t\"params\": [1, 2, -1, 0, 0, 0, 0],\n\t\"type\": \"SimpleItem\"\n}, ";
    // string str5 = "\t\"frame\": 2,\n\t\"params\": [1, 2, -1, 0, 0, 0, 0],\n\t\"type\": \"SimpleItem\"\n} ";
    // string str6 = "\t\t],\n\t\t\"plannedHomePosition\": [30.7461382,103.9291375,519.975],\n\t\t\"vehicleType\": 2,\n\t\t\"version\": 2\n\t},\n\t\"rallyPoints\": {\n\t\t\"points\": [\n\t\t],\n\t\t\"version\": 2\n\t},\n\t\"version\": 1\n}";
    // plan << str0;
    
    // for(int i = 1; i < data_x.size(); i++)
    // {
    //     // cout << data_x[i] << "  " << data_y[i] << endl;
    //     int width = 4;
        
    //     plan << str1 << 2*i << "," << endl;
    //     plan << str2 << data_x[i]/DEG2RAD << ", " << data_y[i]/DEG2RAD << ", " << 15 << "]," << endl;
    //     plan << "\t\"type\": \"SimpleItem\"" << endl << " }," << endl;
    //     plan << str3 << 2*i + 1 << "," << endl;
    //     if(i == (data_x.size() - 1))
    //       plan << str5 << endl;
    //     else
    //     {
    //       plan << str4 << endl;
    //     }
        
    //     // double x = (data_x[i] - origin_lat) * C_EARTH ;
    //     // double y = (data_y[i] - origin_lon) * C_EARTH * cos(data_x[i]) ; 
    //     // cout << "x: " << x << "\t" << "y: " << y << endl;
    // }
    // plan << str6;

    // plan.close();


  //No.2 双扭线形‘8’字
  // int kNumObservations = 100; //100对点   
  // double data[2*kNumObservations];
  // double theta;
  // double ro;
  // for(int i = 0; i < kNumObservations; i++)
    // {
    //     theta = -pi/4 + 0.01*i + 0.01;
    //     ro = sqrt(cos(2*theta));
    //     data[2*i] = a*ro*cos(theta) + c;
    //     data[2*i+1] = b*ro*sin(theta) + d;
    // }
    PLAN test;
    double alt , lat , lon;
    alt = lat = lon = 1;
    test.creatHead();
    test.setHomePosition(origin_lat,origin_lon,530);
    for(int i = 0; i < data_x.size(); i++)
    {
      test.addWayPoint(i + 1 , data_x[i]/DEG2RAD, data_y[i]/DEG2RAD, center_alt);
    }

    // test.printJson();
    test.outJson("waypoint.plan");
    return 0;

}