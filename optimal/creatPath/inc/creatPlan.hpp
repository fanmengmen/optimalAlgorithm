#pragma once

/*建立一个自动生成.plan文件的程序
plan文件属于json文件类型
我使用hlohmann的c++json库，用接近STL的方法实现plan文件的生成*/

#include<iostream>
#include<string>
#include<nlohmann/json.hpp>
#include <iomanip>
#include <fstream>
using json = nlohmann::json;
using namespace std;


class PLAN
{
    public:
        /*************
         ***添加航点***
         *************/
        bool addWayPoint(int num, double lat, double lon, double alt);
        /*设置返航点位置*/
        int setHomePosition(const double home_lat,const double home_lon,const double home_alt);
        /*默认版本12，悬停速度3 巡航速度15 */
        int creatHead(const double hoverSpeed = 3, const double cruiseSpeed = 15, const double firmwareType = 12);
        // 打印json文件
        void printJson();
        // 输出json文件
        int outJson(const string fileName);
    private:
        string str;
        json plan;
};

