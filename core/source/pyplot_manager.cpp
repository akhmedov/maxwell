//
//  pyplot_manager.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 15.10.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "pyplot_manager.hpp"

const std::string PyPlotManager::MULTILINE = 
"import matplotlib.pyplot as plot\n"
"plot.title('$TITLE', fontsize=$FONT)\n"
"plot.xlabel('$XLABEL', fontsize=$FONT)\n"
"plot.ylabel('$YLABEL', fontsize=$FONT)\n"
"$DATA\n"
"plot.legend()\n"
"plot.grid(color ='lightgrey', linestyle='dashed', linewidth=1)\n"
"plot.savefig('$SCRIPT', dpi=300, format='png', bbox_inches='tight')\n"
;

const std::string PyPlotManager::SURFACE =
"import numpy as np\n"
"from matplotlib import cm\n"
"import matplotlib.pyplot as plot\n"
"from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import\n"
"# from matplotlib.ticker import LinearLocator, FormatStrFormatter\n"
"fig = plot.figure()\n"
"ax = fig.gca(projection='3d')\n"
"plot.xlabel('$XLABEL', fontsize=$FONT)\n"
"plot.ylabel('$YLABEL', fontsize=$FONT)\n"
"plot.title('$TITLE', fontsize=$FONT)\n"
"Y, X = np.mgrid[$YFROM:$YTO+$YSTEP:$YSTEP, $XFROM:$XTO+$XSTEP:$XSTEP]\n"
"Z = np.array([\n"
"$DATA\n"
"], dtype=np.float32)\n"
"surf = ax.plot_surface(X, Y, Z, cmap=cm.$COLORMAP, linewidth=0, antialiased=False)\n"
"# ax.set_zlim(-1.01, 1.01)\n"
"# ax.zaxis.set_major_locator(LinearLocator(10))\n"
"# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))\n"
"fig.colorbar(surf, shrink=0.5, aspect=5)\n"
"plot.show()\n"
;

const std::string PyPlotManager::CMAP =
"import numpy as np\n"
"from matplotlib import cm\n"
"from matplotlib import axis as axis\n"
"from matplotlib import pyplot as plot\n"
"from matplotlib import colors as colors\n"
"data = np.array([\n"
"$DATA\n"
"], dtype=np.float32)\n"
"Y, X = np.mgrid[$YFROM:$YTO+$YSTEP:$YSTEP, $XFROM:$XTO+$XSTEP:$XSTEP]\n"
"plot.pcolormesh(X, Y, data, cmap=cm.$COLORMAP)\n"
"plot.xlabel('$XLABEL', fontsize=$FONT)\n"
"plot.ylabel('$YLABEL', fontsize=$FONT)\n"
"plot.title('$TITLE', fontsize=$FONT)\n"
"plot.colorbar()\n"
"plot.savefig('$SCRIPT', dpi=300, format='png', bbox_inches='tight')\n"
;

PyPlotManager::PyPlotManager (const std::string &filename, Logger* global_log)
: ScriptManager(filename, global_log) { }

void PyPlotManager::plot2d (const std::vector<double> &arg, const std::vector<double> &fnc) const
{
	std::string PYPLOT = "plot.plot(X, Y, color='$LCOLOR', linewidth=1, linestyle='$LSTYLE', label='$LLABEL')";
	PYPLOT = std::regex_replace(PYPLOT, std::regex("\\$LLABEL"), this->title);
	PYPLOT = std::regex_replace(PYPLOT, std::regex("\\$LSTYLE"), "solid");

	if (arg.size() != fnc.size()) throw std::invalid_argument("Invalid data: size of arg and fnc must matches.");
	if (arg.empty()) throw std::invalid_argument("Empty plot dataset.");

	std::string X = "X = [", Y = "Y = [";
	for (size_t i = 0; i < fnc.size(); i++) {
		X += (std::isnan(arg[i]) ? "float('nan')" : std::to_string(arg[i])) + ", ";
		Y += (std::isnan(fnc[i]) ? "float('nan')" : std::to_string(fnc[i])) + ", ";
	}
	X += "]\n"; Y += "]\n";

	std::string text = PyPlotManager::MULTILINE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), X+Y+PYPLOT);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy);
	text = std::regex_replace(text, std::regex("\\$TITLE"), this->title);
	text = std::regex_replace(text, std::regex("\\$SCRIPT"), this->fname + ".png");

	switch (this->color) {
		case Colormap::grey: 
			text = std::regex_replace(text, std::regex("\\$LCOLOR"), "black");
			break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	PyPlotManager::write_script(text);
}

void PyPlotManager::plot2d (const std::vector<std::vector<double>> &arg, const std::vector<std::vector<double>> &fnc, const std::vector<std::string> &titles) const
{
	if (arg.size() != fnc.size()) throw std::invalid_argument("Invalid data: size of arg, fnc and titles must matches.");
	if (titles.size() != fnc.size()) throw std::invalid_argument("Invalid data: size of arg, fnc and titles must matches.");
	if (arg.empty()) throw std::invalid_argument("Empty plot dataset.");

	std::string data;

	for (size_t graph = 0; graph < arg.size(); graph++) {

		std::string Xi = "X" + std::to_string(graph);
		std::string Yi = "Y" + std::to_string(graph);
		std::string X = Xi + " = [";
		std::string Y = Yi + " = [";
		for (size_t i = 0; i < fnc[graph].size(); i++) {
			X += (std::isnan(arg[graph][i]) ? "float('nan')" : std::to_string(arg[graph][i])) + ", ";
			Y += (std::isnan(fnc[graph][i]) ? "float('nan')" : std::to_string(fnc[graph][i])) + ", ";
		}
		X += "]\n"; Y += "]\n";

		std::string PYPLOT = "plot.plot($Xi, $Yi, color='$LCOLOR', linewidth=1, linestyle='$LSTYLE', label='$LLABEL')\n";
		PYPLOT = std::regex_replace(PYPLOT, std::regex("\\$Xi"), Xi);
		PYPLOT = std::regex_replace(PYPLOT, std::regex("\\$Yi"), Yi);
		PYPLOT = std::regex_replace(PYPLOT, std::regex("\\$LLABEL"), titles[graph]);
		PYPLOT = std::regex_replace(PYPLOT, std::regex("\\$LSTYLE"), "solid");

		data += X + Y + PYPLOT;
	}
	data = data.substr(0, data.size()-1);

	std::string text = PyPlotManager::MULTILINE;
	text = std::regex_replace(text, std::regex("\\$FONT"), this->font);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy);
	text = std::regex_replace(text, std::regex("\\$TITLE"), this->title);
	text = std::regex_replace(text, std::regex("\\$SCRIPT"), this->fname + ".png");

	switch (this->color) {
		case Colormap::grey: 
			text = std::regex_replace(text, std::regex("\\$LCOLOR"), "black");
			break;
		default: 
			throw std::logic_error("This colorshem is not implemented in GnuPlot::plot_colormap()");
	}

	PyPlotManager::write_script(text);
}

void PyPlotManager::plot3d (const std::vector<std::pair<double,double>> & args, const std::vector<double> & fnc) const
{
	if (args.empty()) throw std::invalid_argument("Empty plot dataset.");
	
	double xmin, xstep, xmax, ymin, ystep, ymax;
	auto matrix = ScriptManager::datagrid_from(args,fnc, xmin, xstep, xmax, ymin, ystep, ymax);

	std::string data;
	for (auto&& line : matrix) {
		data += "\t[";
		for (auto&& val : line)
			data += (std::isnan(val) ? "float('nan')" : std::to_string(val)) + ", ";
		data += "],\n";
	}
	data = data.substr(0, data.size()-1);

	std::string text = PyPlotManager::SURFACE;

	text = std::regex_replace(text, std::regex("\\$XFROM"), std::to_string(xmin));
	text = std::regex_replace(text, std::regex("\\$XSTEP"), std::to_string(xstep));
	text = std::regex_replace(text, std::regex("\\$XTO"), std::to_string(xmax));
	text = std::regex_replace(text, std::regex("\\$YFROM"), std::to_string(ymin));
	text = std::regex_replace(text, std::regex("\\$YSTEP"), std::to_string(ystep));
	text = std::regex_replace(text, std::regex("\\$YTO"), std::to_string(ymax));

	text = std::regex_replace(text, std::regex("\\$FONT"), "16");
	text = std::regex_replace(text, std::regex("\\$SCRIPT"), this->fname + ".png");
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy);
	text = std::regex_replace(text, std::regex("\\$TITLE"), this->title);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);

	switch (this->color) {
		case ScriptManager::Colormap::grey:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "binary");
			break;
		case ScriptManager::Colormap::jet:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "jet");
			break;
		case ScriptManager::Colormap::hot:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "hot");
			break;
		case ScriptManager::Colormap::coolwarm:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "coolwarm");
			break;
		default: 
			throw std::logic_error("This colorshem is not implemented in PyPlotManager::colormap()");
	}

	PyPlotManager::write_script(text);
}

void PyPlotManager::colormap (const std::vector<std::pair<double,double>> & args, const std::vector<double> & fnc) const
{
	if (args.empty()) throw std::invalid_argument("Empty plot dataset.");
	
	double xmin, xstep, xmax, ymin, ystep, ymax;
	auto matrix = ScriptManager::datagrid_from(args,fnc, xmin, xstep, xmax, ymin, ystep, ymax);

	std::string data;
	for (auto&& line : matrix) {
		data += "\t[";
		for (auto&& val : line)
			data += (std::isnan(val) ? "float('nan')" : std::to_string(val)) + ", ";
		data += "],\n";
	}
	data = data.substr(0, data.size()-1);

	std::string text = PyPlotManager::CMAP;

	text = std::regex_replace(text, std::regex("\\$XFROM"), std::to_string(xmin));
	text = std::regex_replace(text, std::regex("\\$XSTEP"), std::to_string(xstep));
	text = std::regex_replace(text, std::regex("\\$XTO"), std::to_string(xmax));
	text = std::regex_replace(text, std::regex("\\$YFROM"), std::to_string(ymin));
	text = std::regex_replace(text, std::regex("\\$YSTEP"), std::to_string(ystep));
	text = std::regex_replace(text, std::regex("\\$YTO"), std::to_string(ymax));

	text = std::regex_replace(text, std::regex("\\$FONT"), "16");
	text = std::regex_replace(text, std::regex("\\$SCRIPT"), this->fname + ".png");
	text = std::regex_replace(text, std::regex("\\$XLABEL"), this->ox);
	text = std::regex_replace(text, std::regex("\\$YLABEL"), this->oy);
	text = std::regex_replace(text, std::regex("\\$TITLE"), this->title);
	text = std::regex_replace(text, std::regex("\\$DATA"), data);

	switch (this->color) {
		case ScriptManager::Colormap::grey:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "binary");
			break;
		case ScriptManager::Colormap::jet:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "jet");
			break;
		case ScriptManager::Colormap::hot:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "hot");
			break;
		case ScriptManager::Colormap::coolwarm:
			text = std::regex_replace(text, std::regex("\\$COLORMAP"), "coolwarm");
			break;
		default: 
			throw std::logic_error("This colorshem is not implemented in PyPlotManager::plot_colormap()");
	}

	PyPlotManager::write_script(text);
}
