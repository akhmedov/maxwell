//
//  gnu_plot.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 18.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//


#include "gnu_plot.hpp"

const std::string GnuPlot::LINE =
"set term x11 font '$FONT'\n"
"set xlabel '$XLABEL'\n"
"set ylabel '$YLABEL'\n"
"set zlabel 'oZ' rotate center\n"
"set grid layerdefault\n"
"$grid << EOD\n"
"$DATA\n"
"EOD\n"
"plot '$grid' using 1:2 with lines linecolor rgb 'black' notitle\n"
"pause -1 'Hit return to continue'\n"
;

const std::string GnuPlot::MULTILINE = 
"set term x11 font '$FONT'\n"
"set xlabel '$XLABEL'\n"
"set ylabel ''\n"
"set zlabel '' rotate center\n"
"set grid layerdefault\n"
"$grid << EOD\n"
"$DATA\n"
"EOD\n"
"$GRAY_LINESTYLE\n"
"plot \\\n"
"$PLOT_LINES\n"
"pause -1 'Hit return to continue'\n"
;

const std::string GnuPlot::SURFACE =
"set term x11 font '$FONT'\n"
"set xlabel '$XLABEL'\n"
"set ylabel '$YLABEL'\n"
"set zlabel '$YLABEL' rotate center\n"
"set grid layerdefault\n"
"$grid << EOD\n"
"$DATA\n"
"EOD\n"
"set hidden3d\n"
"set dgrid3d 50, 50, 50\n"
"$GRAY set style data lines\n"
"$GRAY splot '$grid' using 1:2:3 with lines linecolor rgb 'black' notitle\n"
"$PALURA set pm3d depthorder border linetype -1 linewidth 0.5\n"
"$PALURA set style fill transparent solid 0.65 border\n"
"$PALURA set palette rgb 21,22,23\n"
"$PALURA splot '$grid' using 1:2:3 with pm3d notitle\n"
"pause -1 'Hit return to continue'\n"
;

const std::string GnuPlot::CMAP =
"set terminal png size 1200,1200\n"
"set output '$SCRIPT.png'\n"
"set border linewidth 0\n"
"set palette grey\n"
"$grid << EOD\n"
"$DATA\n"
"EOD\n"
"unset tics\n"
"unset colorbox\n"
"set size square\n"
"$COLOR1 set palette defined (0 1 1 1, 1 0 0 0)\n"
"$COLOR2 set palette defined (0 0 0 0, 1 1 1 0, 2 1 0 0)\n"
"set pm3d map\n"
"set pm3d interpolate 10,10\n"
"splot '$grid' matrix\n"
"# plot '$grid' matrix with image\n"
"pause -1 'Hit return to continue'\n"
;

GnuPlot::GnuPlot (std::string filename)
{
	this->script_name = filename;	

	std::string extention = filename.substr(filename.find_last_of(".") + 1);
	if ( (extention.compare("gnp") != 0) && (extention.compare("txt") != 0))
		this->script_name.append(".gnp");

	this->set_ox_label("oX");
	this->set_oy_label("oY");
	this->set_oz_label("oZ");

	this->font = "times-roman,FONT_SIZE,normal";
	this->set_font_size(20);

	this->is_title_on = false;
	this->is_ox_logscale = false;
	this->is_oy_logscale = false;
	this->is_grid_on = false;
	this->is_plot_in_cage = false;
	this->is_3d_plot = false;
	this->is_colormap_plot = false;

	this->color_schem = Colormap::gray;
	this->gnuplot_path = "gnuplot";
}

void GnuPlot::set_gnuplot_bin (std::string filename)
{
	this->gnuplot_path = filename;
}

void GnuPlot::set_font_size (std::size_t size)
{
	font.replace(
		this->font.find("FONT_SIZE"), 
		std::string("FONT_SIZE").length(), 
		std::to_string(size)
	);
}

void GnuPlot::set_ox_label (std::string text)
{
	this->ox_label = text;
}

void GnuPlot::set_oy_label (std::string text)
{
	this->oy_label = text;
}

void GnuPlot::set_oz_label (std::string text)
{
	this->oz_label = text;
}

void GnuPlot::set_title (std::string text)
{
	this->title = text;
	this->is_title_on = true;
}

void GnuPlot::set_logscale_ox (bool status)
{
	this->is_ox_logscale = status;
}

void GnuPlot::set_logscale_oy (bool status)
{
	this->is_oy_logscale = status;
}

void GnuPlot::grid_on (bool status)
{
	this->is_grid_on = status;
}

void GnuPlot::cage_on ( bool status)
{
	this->is_plot_in_cage = status;
}

void GnuPlot::write_script (const std::string &text) const
{
	std::cout << "Write comands to " << this->script_name << " script... ";
	std::ofstream script;
	script.open( this->script_name );
	script << text;
	script.close();
	std::cout << "Done." << std::endl;
}

std::vector<std::vector<double>> GnuPlot::datagrid_from (const std::vector<std::vector<double>>& datalist, int axis1, int axis2)
{
	std::map<double, std::map<double,double>> maped; // y -> {x -> F}

	for (auto i : datalist) {
		double x = i[axis1];
		double y = i[axis2];
		double F = i.back();
		maped[y][x] = F;
	}

	std::vector<std::vector<double>> grid;

	while (!maped.empty()) {
		
		auto sallestY = maped.begin();
		std::map<double,double> map_sameY = sallestY->second; // get submap with same Y
		maped.erase(sallestY->first); // erase submap from master map

		std::vector<double> vec_sameY; // transform map to vector
		for (auto it = map_sameY.begin(); it != map_sameY.end(); ++it)
        	vec_sameY.push_back( it->second );

		grid.push_back(vec_sameY); // add line to result matrix
	}

	return grid;
}

void GnuPlot::plot_colormap (const std::vector<std::vector<double>> &array, int axis1, int axis2)
{
	auto matrix = GnuPlot::datagrid_from(array,axis1,axis2);

	std::string data;
	for (auto&& line : matrix) {
		for (auto&& val : line)
			data += std::to_string(val) + " ";
		data += "\n";
	}

	std::string text = GnuPlot::CMAP;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$SCRIPT"), this->script_name);

	switch (this->color_schem) {
		case Colormap::gray: {
			text = std::regex_replace(text, std::regex("\\$COLOR2"), "#");
			text = std::regex_replace(text, std::regex("\\$COLOR1"), "");
			} break;
		case Colormap::parula: {
			text = std::regex_replace(text, std::regex("\\$COLOR1"), "#");
			text = std::regex_replace(text, std::regex("\\$COLOR2"), "");
			} break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	this->write_script(text);
}

void GnuPlot::plot2d (const std::vector<std::vector<double>> &array) 
{
	if (array.empty()) throw std::invalid_argument("Empty plot dataset.");
	this->is_3d_plot = false;
	this->is_multi_plot = false;

	std::string data;
	for (auto&& line : array) {
		for (auto&& val : line)
			data += std::to_string(val) + " ";
		data += "\n";
	}

	std::string text = GnuPlot::LINE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox_label);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy_label);

	switch (this->color_schem) {
		case Colormap::gray: { } break;
		case Colormap::parula: { } break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	this->write_script(text);
}

void GnuPlot::plot_multi (const std::vector<std::vector<double>> &arrays, const std::vector<std::string> &title) 
{
	std::string LT_FOR = "set for [i=1:$LINES] linetype i dt i\n";
	const std::string LT_ITEM = "set style line $ITEM lt $ITEM lc rgb 'black' lw 1 \n";
	const std::string PLOT = "'$grid' using 1:1+$ITEM with lines ls $ITEM title '$TITLE',\\\n";

	if (arrays.empty()) throw std::invalid_argument("Empty plot dataset.");
	if (arrays[0].size()-1 != title.size()) throw std::invalid_argument("Size of input arrays does not match.");

	std::string data;
	for (auto&& line : arrays) {
		for (auto&& val : line)
			data += std::to_string(val) + " ";
		data += "\n";
	}

	std::string text = GnuPlot::MULTILINE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox_label);

	switch (this->color_schem) {
		case Colormap::gray: {
			std::string linetype = std::regex_replace(LT_FOR, std::regex("\\$LINES"), std::to_string(title.size()));
			for (std::size_t i = 1; i <= title.size(); i++) {
				std::string lti = LT_ITEM;
				lti = std::regex_replace(lti, std::regex("\\$ITEM"), std::to_string(i));
				linetype.append(lti);
			}
			text = std::regex_replace(text, std::regex("\\$GRAY_LINESTYLE"), linetype);
		} break;
		case Colormap::parula: { } break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	std::string plot_cmd;
	for (std::size_t i = 1; i <= title.size(); i++) {
		std::string plot = PLOT;
		plot = std::regex_replace(plot, std::regex("\\$ITEM"), std::to_string(i));
		plot = std::regex_replace(plot, std::regex("\\$TITLE"), title[i-1]);
		plot_cmd.append(plot);
	}

	plot_cmd = plot_cmd.substr(0, plot_cmd.size()-3);
	text = std::regex_replace(text, std::regex("\\$PLOT_LINES"), plot_cmd);
	
	this->write_script(text);
}

void GnuPlot::plot3d (const std::vector<std::vector<double>> &matrix)
{
	if (matrix.empty()) throw std::invalid_argument("Empty plot dataset.");

	std::string data;
	for (auto&& point : matrix) {
		std::string x = std::to_string(point[1]).append(" ");
		std::string y = std::to_string(point[2]).append(" ");
		std::string z = std::to_string(point[0]).append(" ");
		std::string l = x.append(y).append(z).append("\n");
		data.append(l);
	}

	std::string text = GnuPlot::SURFACE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox_label);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy_label);
	text = std::regex_replace(text, std::regex("\\$ZLABEL"), this->oz_label);

	switch (this->color_schem) {
		case Colormap::gray: {
			text = std::regex_replace(text, std::regex("\\$PALURA"), "#");
			text = std::regex_replace(text, std::regex("\\$GRAY"), "");
		} break;
		case Colormap::parula: {
			text = std::regex_replace(text, std::regex("\\$PALURA"), "");
			text = std::regex_replace(text, std::regex("\\$GRAY"), "#");
		} break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	this->write_script(text);
}

void GnuPlot::call_gnuplot ()
{
	std::cout << "Call gnuplot binary..." << std::endl;
	std::cout.flush();
	std::string command = this->gnuplot_path.append(" -c ");
	system( command.append(this->script_name).c_str() );
}

void GnuPlot::set_colormap (const Colormap &schem)
{
	color_schem = schem;
}

/* void GnuPlot::direct_gnuplot_call (const Text &plot_data) const
{
	std::string cmd_list;
	cmd_list.append("set term x11; ");
	// if (this->is_title_on) cmd_list.("set title \"" << this->title << "\";");
	if (this->is_3d_plot) thi
} */
