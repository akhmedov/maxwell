//
//  gnuplot_manager.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 18.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "gnuplot_manager.hpp"

const std::string GnuPlotManager::LINE =
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

const std::string GnuPlotManager::MULTILINE = 
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

const std::string GnuPlotManager::SURFACE =
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

const std::string GnuPlotManager::CMAP =
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

GnuPlotManager::GnuPlotManager (const std::string &filename, Logger* global_log)
: ScriptManager(filename, global_log) { }

void GnuPlotManager::plot2d (const std::vector<double> &arg, const std::vector<double> &fnc) const
{
	if (arg.size() != fnc.size()) throw std::invalid_argument("Invalid data: size of arg and fnc must matches.");
	if (arg.empty()) throw std::invalid_argument("Empty plot dataset.");

	std::string data;
	for (size_t i = 0; i < fnc.size(); i++)
		data += std::to_string(arg[i]) + " " + std::to_string(fnc[i]) + "\n";

	std::string text = GnuPlotManager::LINE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy);

	switch (this->color) {
		case Colormap::grey: { } break;
		case Colormap::hot: { } break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlotManager::plot_colormap()");
	}

	this->write_script(text);
}

void GnuPlotManager::plot2d (const std::vector<std::vector<double>> &arg, const std::vector<std::vector<double>> &fnc, const std::vector<std::string> &titles) const
{
	if (arg.size() != fnc.size()) throw std::invalid_argument("Invalid data: size of arg, fnc and titles must matches.");
	if (titles.size() != fnc.size()) throw std::invalid_argument("Invalid data: size of arg, fnc and titles must matches.");
	if (arg.empty()) throw std::invalid_argument("Empty plot dataset.");

	std::string LT_FOR = "set for [i=1:$LINES] linetype i dt i\n";
	const std::string LT_ITEM = "set style line $ITEM lt $ITEM lc rgb 'black' lw 1 \n";
	const std::string PLOT = "'$grid' using 1:1+$ITEM with lines ls $ITEM title '$TITLE',\\\n";

	// BUG: only arg[0] array is used for gnuplot - different argument is not suported
	std::string data;
	for (size_t vert = 0; vert < arg[0].size(); vert++) {
		data += std::to_string(arg[0][vert]) + " ";
		for (size_t horiz = 0; horiz < fnc.size(); horiz++) {
			data += std::to_string(fnc[horiz][vert]) + " ";
		}
		data += "\n";
	}

	std::string text = GnuPlotManager::MULTILINE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox);

	switch (this->color) {
		case Colormap::grey: {
			std::string linetype = std::regex_replace(LT_FOR, std::regex("\\$LINES"), std::to_string(titles.size()));
			for (std::size_t i = 1; i <= titles.size(); i++) {
				std::string lti = LT_ITEM;
				lti = std::regex_replace(lti, std::regex("\\$ITEM"), std::to_string(i));
				linetype.append(lti);
			}
			text = std::regex_replace(text, std::regex("\\$GRAY_LINESTYLE"), linetype);
		} break;
		case Colormap::hot: { } break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlotManager::plot_colormap()");
	}

	std::string plot_cmd;
	for (std::size_t i = 1; i <= titles.size(); i++) {
		std::string plot = PLOT;
		plot = std::regex_replace(plot, std::regex("\\$ITEM"), std::to_string(i));
		plot = std::regex_replace(plot, std::regex("\\$TITLE"), titles[i-1]);
		plot_cmd.append(plot);
	}

	plot_cmd = plot_cmd.substr(0, plot_cmd.size()-3);
	text = std::regex_replace(text, std::regex("\\$PLOT_LINES"), plot_cmd);
	
	this->write_script(text);
}

void GnuPlotManager::plot3d (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const
{
	if (args.empty()) throw std::invalid_argument("Empty plot dataset.");
	
	double xmin, xstep, xmax, ymin, ystep, ymax; // unused
	auto matrix = ScriptManager::datagrid_from(args,fnc, xmin, xstep, xmax, ymin, ystep, ymax);

	std::string data;
	for (auto&& point : matrix) {
		std::string x = std::to_string(point[1]).append(" ");
		std::string y = std::to_string(point[2]).append(" ");
		std::string z = std::to_string(point[0]).append(" ");
		std::string l = x.append(y).append(z).append("\n");
		data.append(l);
	}

	std::string text = GnuPlotManager::SURFACE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy);
	text = std::regex_replace(text, std::regex("\\$ZLABEL"), "// TODO:");

	switch (this->color) {
		case Colormap::grey: {
			text = std::regex_replace(text, std::regex("\\$PALURA"), "#");
			text = std::regex_replace(text, std::regex("\\$GRAY"), "");
		} break;
		case Colormap::hot: {
			text = std::regex_replace(text, std::regex("\\$PALURA"), "");
			text = std::regex_replace(text, std::regex("\\$GRAY"), "#");
		} break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	this->write_script(text);
}

void GnuPlotManager::colormap (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const
{
	if (args.empty()) throw std::invalid_argument("Empty plot dataset.");
	
	double xmin, xstep, xmax, ymin, ystep, ymax; // unused
	auto matrix = ScriptManager::datagrid_from(args,fnc, xmin, xstep, xmax, ymin, ystep, ymax);

	std::string data;
	for (auto&& line : matrix) {
		for (auto&& val : line)
			data += std::to_string(val) + " ";
		data += "\n";
	}

	std::string text = GnuPlotManager::CMAP;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$SCRIPT"), this->fname);

	switch (this->color) {
		case Colormap::grey: {
			text = std::regex_replace(text, std::regex("\\$COLOR2"), "#");
			text = std::regex_replace(text, std::regex("\\$COLOR1"), "");
			} break;
		case Colormap::hot: {
			text = std::regex_replace(text, std::regex("\\$COLOR1"), "#");
			text = std::regex_replace(text, std::regex("\\$COLOR2"), "");
			} break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	this->write_script(text);
}
