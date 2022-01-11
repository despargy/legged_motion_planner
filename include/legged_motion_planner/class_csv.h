#ifndef __CLASS_CSV_H__
#define __CLASS_CSV_H__
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen/Dense>

using namespace std;

struct labels{

  Eigen::Vector3d CoM_positions;
  Eigen::Vector3d CoM_linear_velocities;
  Eigen::Vector3d CoM_linear_accelerations;

  Eigen::Vector3d footL_positions;
  Eigen::Vector3d footR_positions;

  Eigen::Vector3d ZMP_positions;

};

vector<labels> pass_csv(string fname)
{

	vector<vector<string>> content;
	vector<string> row;
	string line, word;

	fstream file (fname, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();

			stringstream str(line);

			while(getline(str, word, ','))
				row.push_back(word);

			content.push_back(row);
		}
	}
	else
		cout<<"Could not open the file\n";


  vector<labels> results;
  struct labels intp;

  for(int i=0;i<content.size();i++)
	{
    int j = 0;

    intp.CoM_positions(0)=stod(content[i][j]); j++;
    intp.CoM_positions(1)=stod(content[i][j]); j++;
     intp.CoM_positions(2)=stod(content[i][j]); j++;

    intp.CoM_linear_velocities(0)=stod(content[i][j]); j++;intp.CoM_linear_velocities(1)=stod(content[i][j]); j++; intp.CoM_linear_velocities(2)=stod(content[i][j]); j++;
    intp.CoM_linear_accelerations(0)=stod(content[i][j]); j++;intp.CoM_linear_accelerations(1)=stod(content[i][j]); j++; intp.CoM_linear_accelerations(2)=stod(content[i][j]); j++;

    intp.footL_positions(0)=stod(content[i][j]); j++;intp.footL_positions(1)=stod(content[i][j]); j++; intp.footL_positions(2)=stod(content[i][j]); j++;
    intp.footR_positions(0)=stod(content[i][j]); j++;intp.footR_positions(1)=stod(content[i][j]); j++; intp.footR_positions(2)=stod(content[i][j]); j++;

    intp.ZMP_positions(0)=stod(content[i][j]); j++;intp.ZMP_positions(1)=stod(content[i][j]); j++; intp.ZMP_positions(2)=stod(content[i][j]);

    results.push_back(intp);
	}

  // for(int i=0;i<results.size();i++)
  // {
  // 	// for(int j=0;j<content[i].size();j++)
  // 	// {
  // 		cout<<results[i].CoM_positions<<" ";
  // 	// }
  // 	cout<<"\n";
  // }

	return results;
}
#endif
