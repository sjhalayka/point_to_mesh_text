#ifndef main_H
#define main_H



#include "marching_cubes.h"




#include <cstdlib>

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;
using std::istringstream;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <set>
using std::set;

#include <map>
using std::map;

#include <utility>
using std::pair;

#include <ios>
using std::ios_base;







bool write_triangles_to_binary_stereo_lithography_file(const vector<triangle>& triangles, const char* const file_name)
{
	cout << "Triangle count: " << triangles.size() << endl;

	if (0 == triangles.size())
		return false;

	// Write to file.
	ofstream out(file_name, ios_base::binary);

	if (out.fail())
		return false;

	const size_t header_size = 80;
	vector<char> buffer(header_size, 0);
	const unsigned int num_triangles = static_cast<unsigned int>(triangles.size()); // Must be 4-byte unsigned int.
	vertex_3 normal;

	// Write blank header.
	out.write(reinterpret_cast<const char*>(&(buffer[0])), header_size);

	// Write number of triangles.
	out.write(reinterpret_cast<const char*>(&num_triangles), sizeof(unsigned int));

	// Copy everything to a single buffer.
	// We do this here because calling ofstream::write() only once PER MESH is going to 
	// send the data to disk faster than if we were to instead call ofstream::write()
	// thirteen times PER TRIANGLE.
	// Of course, the trade-off is that we are using 2x the RAM than what's absolutely required,
	// but the trade-off is often very much worth it (especially so for meshes with millions of triangles).
	cout << "Generating normal/vertex/attribute buffer" << endl;

	// Enough bytes for twelve 4-byte floats plus one 2-byte integer, per triangle.
	const size_t data_size = (12 * sizeof(float) + sizeof(short unsigned int)) * num_triangles;
	buffer.resize(data_size, 0);

	// Use a pointer to assist with the copying.
	// Should probably use std::copy() instead, but memcpy() does the trick, so whatever...
	char* cp = &buffer[0];

	for (vector<triangle>::const_iterator i = triangles.begin(); i != triangles.end(); i++)
	{
		// Get face normal.
		vertex_3 v0 = i->vertex[1] - i->vertex[0];
		vertex_3 v1 = i->vertex[2] - i->vertex[0];

		normal = v0.cross(v1);
		normal.normalize();
		normal.x = -normal.x;
		normal.y = -normal.y;
		normal.z = -normal.z;

		memcpy(cp, &normal.x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &normal.y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &normal.z, sizeof(float)); cp += sizeof(float);

		memcpy(cp, &i->vertex[2].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[2].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[2].z, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[1].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[1].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[1].z, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[0].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[0].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[0].z, sizeof(float)); cp += sizeof(float);

		cp += sizeof(short unsigned int);
	}

	cout << "Writing " << data_size / 1048576.0f << " MB of data to binary Stereo Lithography file: " << file_name << endl;

	out.write(reinterpret_cast<const char*>(&buffer[0]), data_size);
	out.close();

	return true;
}



void tesselate_field(const vector<float>& values, vector<triangle>& triangle_list, const float& isovalue, const float& grid_min, const float& grid_max, const size_t& res)
{

	triangle_list.clear();

	const float step_size = (grid_max - grid_min) / (res - 1);

	for (size_t x = 0; x < res - 1; x++)
	{
		for (size_t y = 0; y < res - 1; y++)
		{
			for (size_t z = 0; z < res - 1; z++)
			{
				grid_cube temp_cube;

				size_t x_offset = 0;
				size_t y_offset = 0;
				size_t z_offset = 0;

				// Setup vertex 0
				x_offset = 0;
				y_offset = 0;
				z_offset = 0;
				temp_cube.vertex[0].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[0].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[0].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[0] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Setup vertex 1
				x_offset = 1;
				y_offset = 0;
				z_offset = 0;
				temp_cube.vertex[1].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[1].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[1].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[1] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Setup vertex 2
				x_offset = 1;
				y_offset = 0;
				z_offset = 1;
				temp_cube.vertex[2].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[2].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[2].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[2] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Setup vertex 3
				x_offset = 0;
				y_offset = 0;
				z_offset = 1;
				temp_cube.vertex[3].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[3].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[3].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[3] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Setup vertex 4
				x_offset = 0;
				y_offset = 1;
				z_offset = 0;
				temp_cube.vertex[4].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[4].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[4].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[4] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Setup vertex 5
				x_offset = 1;
				y_offset = 1;
				z_offset = 0;
				temp_cube.vertex[5].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[5].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[5].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[5] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Setup vertex 6
				x_offset = 1;
				y_offset = 1;
				z_offset = 1;
				temp_cube.vertex[6].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[6].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[6].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[6] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Setup vertex 7
				x_offset = 0;
				y_offset = 1;
				z_offset = 1;
				temp_cube.vertex[7].x = grid_min + ((x + x_offset) * step_size);
				temp_cube.vertex[7].y = grid_min + ((y + y_offset) * step_size);
				temp_cube.vertex[7].z = grid_min + ((z + z_offset) * step_size);
				temp_cube.value[7] = values[(x + x_offset) * (res) * (res)+(y + y_offset) * (res)+(z + z_offset)];

				// Generate triangles from cube
				triangle temp_triangle_array[5];

				short unsigned int number_of_triangles_generated = tesselate_grid_cube(isovalue, temp_cube, temp_triangle_array);

				for (short unsigned int i = 0; i < number_of_triangles_generated; i++)
				{
					triangle_list.push_back(temp_triangle_array[i]);
				}
			}
		}
	}
}


	
void blur_field(vector<float> &input_field, size_t res)
{
	vector<float> field(res * res * res, 0.0);

	for (size_t i = 1; i < res - 1; i++)
	{
		for (size_t j = 1; j < res - 1; j++)
		{
			for (size_t k = 1; k < res - 1; k++)
			{
				float sum = 0;

				for (signed char l = -1; l <= 1; l++)
				{
					for (signed char m = -1; m <= 1; m++)
					{
						for (signed char n = -1; n <= 1; n++)
						{
							size_t index = (k + n) * res * res;
							index += (j + m) * res;
							index += (i + l);

							sum += input_field[index];
						}
					}
				}

				size_t centre_index = k * res * res;
				centre_index += j * res;
				centre_index += i;

				field[centre_index] = sum / 27.0f;
			}
		}
	}

	input_field = field;
}



void convert_point_cloud_to_mesh(const char* const points_filename, size_t res, const char* const stl_filename)
{
	vector<float> field;
	field.resize(res * res * res, 0.0);

	ifstream in_file(points_filename);

	float curr_x_min = numeric_limits<float>::max();
	float curr_y_min = numeric_limits<float>::max();
	float curr_z_min = numeric_limits<float>::max();
	float curr_x_max = -numeric_limits<float>::max();
	float curr_y_max = -numeric_limits<float>::max();
	float curr_z_max = -numeric_limits<float>::max();

	string line;

	size_t count = 0;

	while(getline(in_file, line))
	{
		if (line == "")
			continue;

		count++;

		if (count % 100000 == 0)
			cout << count / static_cast<long double>(res * res * res) << endl;

		istringstream iss(line);

		float x, y, z;
		size_t num;

		iss >> x;
		iss >> y;
		iss >> z;
		
		num = 8;

		if (x < curr_x_min)
			curr_x_min = x;

		if (x > curr_x_max)
			curr_x_max = x;

		if (y < curr_y_min)
			curr_y_min = y;

		if (y > curr_y_max)
			curr_y_max = y;

		if (z < curr_z_min)
			curr_z_min = z;

		if (z > curr_z_max)
			curr_z_max = z;
	}

	float x_extent = curr_x_max - curr_x_min;
	float y_extent = curr_y_max - curr_y_min;
	float z_extent = curr_z_max - curr_z_min;

	in_file.close();
	in_file.open(points_filename);

	count = 0;

	while (getline(in_file, line))
	{
		if (line == "")
			continue;

		count++;

		if (count % 100000 == 0)
			cout << count / static_cast<long double>(res * res * res) << endl;

		istringstream iss(line);

		float x, y, z;
		size_t num;

		iss >> x;
		iss >> y;
		iss >> z;
		
		num = 8;

		float x_location = x - curr_x_min;
		size_t x_index = static_cast<size_t>(static_cast<double>(res) * (x_location / x_extent));

		if (x_index >= res)
			x_index = res - 1;

		float y_location = y - curr_y_min;
		size_t y_index = static_cast<size_t>(static_cast<double>(res) * (y_location / y_extent));

		if (y_index >= res)
			y_index = res - 1;

		float z_location = z - curr_z_min;
		size_t z_index = static_cast<size_t>(static_cast<double>(res) * (z_location / z_extent));

		if (z_index >= res)
			z_index = res - 1;

		size_t index = z_index * res * res;
		index += y_index * res;
		index += x_index;

		if (x_index == 0 || x_index == res - 1 ||
			y_index == 0 || y_index == res - 1 ||
			z_index == 0 || z_index == res - 1)
		{
			field[index] = 0; // add blank border
		}
		else
		{
			field[index] += static_cast<float>(num);
		}
	}

	blur_field(field, res);
	blur_field(field, res);
	blur_field(field, res);
	blur_field(field, res);
	blur_field(field, res);


	vector<triangle> triangles;
	tesselate_field(field, triangles, 0.1f, curr_x_min, curr_x_max, res);
	write_triangles_to_binary_stereo_lithography_file(triangles, stl_filename);
}





#endif
