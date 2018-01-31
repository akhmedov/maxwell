#include "gnu_plot.hpp" 

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

void GnuPlot::plot2d (const std::vector<std::vector<double>> &array) 
{
	this->is_3d_plot = false;
	this->is_multi_plot = false;
	
	std::vector<std::string> plot_data;
	for (auto&& point : array) {
		std::string str_line;
		std::string x = std::to_string(point[0]).append(" ");
		std::string y = std::to_string(point[1]).append(" ");
		str_line.append( x.append(y) );
		plot_data.push_back(str_line);
	}
	
	std::cout << "Write comands to " << this->script_name << " script... ";
	std::cout.flush();
	this->write_commants_to_script(plot_data);
	std::cout << "Done." << std::endl;
	std::cout.flush();
}

void GnuPlot::plot_multi (const std::vector<std::vector<double>> &arrays, const std::vector<std::string> &title) 
{
	this->is_3d_plot = false;
	this->is_multi_plot = true;
	
	std::vector<std::string> plot_data;
	for (auto&& line : arrays) {
		std::string str_line;
		for (auto&& value : line) {
			std::string val = std::to_string(value).append(" ");
			str_line.append(val);
		}
		plot_data.push_back(str_line);
	}
	
	std::cout << "Write comands to " << this->script_name << " script... ";
	std::cout.flush();
	this->write_commants_to_script(plot_data, title);
	std::cout << "Done." << std::endl;
	std::cout.flush();
}

void GnuPlot::plot3d (const std::vector<std::vector<double>> &matrix)
{
	this->is_3d_plot = true;
	this->is_multi_plot = false;

	std::vector<std::string> plot_data;
	for (auto&& point : matrix) {
		std::string str_line;
		std::string x = std::to_string(point[1]).append(" ");
		std::string y = std::to_string(point[2]).append(" ");
		std::string z = std::to_string(point[0]).append(" ");
		str_line.append( x.append(y).append(z) );
		plot_data.push_back(str_line);
	}

	std::cout << "Write comands to " << this->script_name << " script... ";
	std::cout.flush();
	this->write_commants_to_script(plot_data);
	std::cout << "Done." << std::endl;
	std::cout.flush();
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

void GnuPlot::write_commants_to_script (const std::vector<std::string> &plot_data, const std::vector<std::string> &title) const
{
	std::ofstream script;
	script.open( this->script_name );
	script << "set term x11 font \"" << this->font << "\" \n";

	if (this->is_title_on) script << "set title \"" << this->title << "\"\n";

	script << "set xlabel \"" << this->ox_label << "\"\n";
	script << "set ylabel \"" << this->oy_label << "\"\n";
	script << "set zlabel \"" << this->oz_label << '\"' << " rotate center" << '\n';

	if (this->is_ox_logscale) script << "set logscale x \n";
	if (this->is_oy_logscale) script << "set logscale y \n";
	if (this->is_grid_on) script << "set grid layerdefault" << '\n';

	script << "$grid << EOD \n";
	for (auto i : plot_data) script << i << '\n';
	script << "EOD \n";

	if (this->is_3d_plot) {
		if (this->color_schem == Colormap::gray) {
			script << "set hidden3d \n"; // pm3d
			script << "set dgrid3d 50, 50, 50 \n";
			script << "set style data lines \n";
			script << "splot '$grid' using 1:2:3 with lines linecolor rgb \"black\" notitle \n";
		} else if (this->color_schem == Colormap::parula) {
			script << "set dgrid3d 50, 50, 50 \n";
			script << "set pm3d depthorder border linetype -1 linewidth 0.5\n";
			script << "set style fill  transparent solid 0.65 border\n";
			script << "set palette rgb 21,22,23\n";
			script << "set hidden3d\n";
			script << "splot '$grid' using 1:2:3 with pm3d notitle \n";
		}

	} else if (this->is_multi_plot) {

		std::size_t curves = std::count_if(plot_data[0].begin(), plot_data[0].end(),
			[] (char c) { return std::isspace(c); }
		) - 1;

		if (this->color_schem == Colormap::gray) {
			script << "set for [i=1:" << std::to_string(curves) << "] linetype i dt i \n";
			for (std::size_t i = 1; i <= curves; i++)
				script << "set style line " << std::to_string(i) << " lt " << std::to_string(i) << " lc rgb \"black\" lw 1 \n";
		}

		std::string plot_conf = "'$grid' using 1:COLUMN";
		std::string options = " with lines ls STYLE title \"LABLE\"";

		for (std::size_t i = 1; i <= curves; i++) {
			std::string prefix = (i == 1) ? "plot " : " ";
			std::string postfix = (i != curves) ? "," : "\n";
			std::string cmd = prefix.append(plot_conf).append(options).append(postfix);
			cmd.replace(cmd.find("STYLE"), std::string("STYLE").length(), std::to_string(i));
			cmd.replace(cmd.find("COLUMN"), std::string("COLUMN").length(), std::to_string(i+1));
			cmd.replace(cmd.find("LABLE"), std::string("LABLE").length(), title[i-1]);
			script << cmd;
		}
	} else {
		script << "plot '$grid' using 1:2 with lines linecolor rgb \"black\" notitle \n";
	}

	script << "pause -1 \"Hit return to continue\" \n";
	script.close();
}
