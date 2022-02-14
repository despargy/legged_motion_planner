#ifndef __CLASS_CSV_H__
#define __CLASS_CSV_H__
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <legged_motion_planner/labels.h>

vector<labels> pass_csv(string fname, int NUM_LEGS)
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
  intp.set_foot_size(NUM_LEGS);

  for(int i=0;i<content.size();i++)
	{
    int j = 0;

    intp.CoM_positions(0)=stod(content[i][j]); j++;
    intp.CoM_positions(1)=stod(content[i][j]); j++;
     intp.CoM_positions(2)=stod(content[i][j]); j++;

    intp.CoM_linear_velocities(0)=stod(content[i][j]); j++;intp.CoM_linear_velocities(1)=stod(content[i][j]); j++; intp.CoM_linear_velocities(2)=stod(content[i][j]); j++;
    intp.CoM_linear_accelerations(0)=stod(content[i][j]); j++;intp.CoM_linear_accelerations(1)=stod(content[i][j]); j++; intp.CoM_linear_accelerations(2)=stod(content[i][j]); j++;

    for(int ll=0; ll < NUM_LEGS; ll++)
    {
      intp.foot_position_v[ll](0)=stod(content[i][j]); j++;
			intp.foot_position_v[ll](1)=stod(content[i][j]); j++;
			intp.foot_position_v[ll](2)=stod(content[i][j]); j++;
    }

    intp.ZMP_positions(0)=stod(content[i][j]); j++;intp.ZMP_positions(1)=stod(content[i][j]); j++; intp.ZMP_positions(2)=stod(content[i][j]);

    results.push_back(intp);
	}

  // for(int i=0;i<results.size();i++)
  // {
  //   cout<<"a result: ";
  //
  // 	for(int ll=0;ll<NUM_LEGS;ll++)
  // 	{
  //     cout<<"a LEG has 3 foot pos: ";
  // 		cout<<results[i].foot_position_v[ll](0)<<" ";
  //     cout<<results[i].foot_position_v[ll](1)<<" ";
  //     cout<<results[i].foot_position_v[ll](2)<<" ";
  //     cout<<"\n";
  // 	}
  // 	cout<<"\n";
  // }

	return results;
}
#endif
