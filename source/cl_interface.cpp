//
//  cl_interface.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 28.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "cl_interface.hpp"

CLI::CLI (Config* conf)
{
	this->global_conf = conf;

	/* TODO: set short_binar_cmd */

	std::vector<std::string> argtype_regs = {
		"INT",
		"FLOAT",
		"PATH",
		"PAIR", /* FLOAT,FLOAT */
		"EVENT" /* FLOAT:FLOAT,FLOAT:FLOAT,FLOAT,FLOAT */
	};

	this->binar_def_cmd = {
		/* impulse shape */
		{"--shape",		"rect"},
		{"--shape",		"gauss"},
		{"--shape",		"triangle"}
	};

	this->binar_var_cmd = {
		/* culculation options */
		{"--plot",		"INT"},
		{"--data",		"INT"},
		{"--conf", 		"PATH"},
		{"--duration", 	"FLOAT"},
		{"--noise",		"FLOAT"},
		/* model options */
		{"--magnitude", "FLOAT"},
		{"--radius", 	"FLOAT"},
		{"--epsr", 		"FLOAT"},
		{"--mur", 		"FLOAT"},
		{"--kerr", 		"FLOAT"},
		/* reciver position */
		{"--cartesian",	"EVENT"},
		{"--cylindric",	"EVENT"}
	};

	this->binar_cmd.insert(binar_def_cmd.begin(),binar_def_cmd.end());
	this->binar_cmd.insert(binar_var_cmd.begin(),binar_var_cmd.end());

	this->unar_cmd = {
		{"--help"},
		{"--version"},
		{"--safe"}
	};
}

void CLI::call_handler (const std::vector<std::string>& argv) const
{
	// Configuration conf = Configuration::Configuration();

	if (argv.size() == 1) {
		std::cout << "Error: invalid usage" << std::endl;
		std::cout << "Call --help for more detales." << std::endl;
		return;
	}

	for (std::size_t iter = 1; iter < argv.size(); iter++) {

		if ( binar_cmd.find(argv[iter]) != binar_cmd.end() ) {

			if (iter == argv.size() - 1) {
				std::cout << "Error: no argument for paramiter ";
				std::cout << argv[iter] << '.' << std::endl;
				std::cout << "Call --help for more detales." << std::endl;
				return;
			}

			bool valid_arg = false;
			auto iter_args = binar_def_cmd.equal_range(argv[iter]);
			for (auto i = iter_args.first; i != iter_args.second; ++i)
				if (!argv[iter+1].compare(i->second)) valid_arg = true;

			iter_args = binar_var_cmd.equal_range(argv[iter]);
			for (auto i = iter_args.first; i != iter_args.second; ++i) {
				if (CLI::is_int   (argv[iter+1]) && !(i->second.compare("INT"  ))) valid_arg = true;
				if (CLI::is_float (argv[iter+1]) && !(i->second.compare("FLOAT"))) valid_arg = true;
				if (CLI::is_path  (argv[iter+1]) && !(i->second.compare("PATH" ))) valid_arg = true;
				if (CLI::is_event (argv[iter+1]) && !(i->second.compare("EVENT"))) valid_arg = true;
			}

			if (!valid_arg)  {
				std::cout << "Error: " << argv[iter+1];
				std::cout << " is not valid argument for paramiter ";
				std::cout << argv[iter] << std::endl;
				std::cout << "Call --help for more detales." << std::endl;
				return;
			}
			
			this->update_config(argv[iter], argv[iter+1]);
			iter++;

		} else if ( unar_cmd.find(argv[iter]) != unar_cmd.end() ) {
			this->update_config(argv[iter]);

		} else {
			std::cout << "Error: " << argv[iter];
			std::cout << " is not valid paramiter or argument." << std::endl;
			std::cout << "Call --help for more detales." << std::endl;
			return;
		}
	}

	this->run_list();
}

void CLI::run_list () const
{
	/* read config file */
	this->read_config_file();

	if (this->global_conf->version_opt()) CLI::print_version();
	else if (this->global_conf->help_opt()) CLI::print_help(); 
	else if (this->global_conf->dataset_model()) {
		throw std::logic_error ("Not implemented!");
		/* CLI::print_arguments();
		data_model->call(global_conf->dataset_model()); */
	} else if (this->global_conf->plot_model()) {
		CLI::print_arguments();
		std::size_t model_name = this->global_conf->plot_model();
		plot_model->call( (PlotModel::Name) model_name);
	}
}

void CLI::update_config (const std::string& param, const std::string& arg) const
{
	if (!param.compare("--shape")) {
		if (!arg.compare("rect")) {
			this->global_conf->impulse_shape(ImpulseShape::rect);
		} else throw std::logic_error("Only rect shape is implemented!");
	}

	if (!param.compare("--plot")) {
		this->global_conf->plot_model(std::stoi(arg));
		FieldComponent comp = (FieldComponent) std::stoi(arg);
		this->global_conf->field_component(comp);
	}

	if (!param.compare("--data"))
		throw std::logic_error("Dataset models are not implemented!");

	if (!param.compare("--conf"))
		this->global_conf->maxwell_config_path(arg);

	if (!param.compare("--duration"))
		throw std::logic_error("Duration option is not implemented!");

	if (!param.compare("--magnitude"))
		this->global_conf->plane_disk_magnitude(std::stod(arg));

	if (!param.compare("--radius"))
		this->global_conf->plane_disk_radius(std::stod(arg));

	if (!param.compare("--epsr"))
		this->global_conf->plane_disk_epsr(std::stod(arg));

	if (!param.compare("--mur"))
		this->global_conf->plane_disk_mur(std::stod(arg));

	if (!param.compare("--kerr")) {
		if (std::stod(arg) != 0) {
			this->global_conf->kerr_medium(true);
			this->global_conf->kerr_value(std::stod(arg));
		}
	}

	if (!param.compare("--noise"))
		this->global_conf->noise_level(std::stod(arg));

	if (!param.compare("--cartesian"))
		// TODO: cartesian convert to cylindric
		throw std::logic_error( "TODO: --cartesian handle is not implemented" );

	if (!param.compare("--cylindric")) {
		std::vector<std::string> wp = CLI::split(arg,',');

		if (CLI::is_range_value(wp[0])) {
			std::vector<std::string> vt = CLI::split(wp[0],':');
			this->global_conf->receiver_vt(std::stod(vt[0]), std::stod(vt[1]), std::stod(vt[2]));
		} else this->global_conf->receiver_vt(std::stod(wp[0]));

		if (CLI::is_range_value(wp[1])) {
			std::vector<std::string> rho = CLI::split(wp[1],':');
			this->global_conf->receiver_rho(std::stod(rho[0]), std::stod(rho[1]), std::stod(rho[2]));
		} else this->global_conf->receiver_rho(std::stod(wp[1]));

		double rad = M_PI / 180;
		if (CLI::is_range_value(wp[2])) {
			std::vector<std::string> phi = CLI::split(wp[2],':');
			double rad0 = rad*std::stod(phi[0]); rad0 = roundf(rad0 * 100000) / 100000;
			double rad1 = rad*std::stod(phi[1]); rad1 = roundf(rad1 * 100000) / 100000;
			double rad2 = rad*std::stod(phi[2]); rad2 = roundf(rad2 * 100000) / 100000;
			this->global_conf->receiver_phi(rad0,rad1,rad2);
		} else this->global_conf->receiver_phi(rad*std::stod(wp[2]));

		if (CLI::is_range_value(wp[3])) {
			std::vector<std::string> z = CLI::split(wp[3],':');
			this->global_conf->receiver_z(std::stod(z[0]), std::stod(z[1]), std::stod(z[2]));
		} else this->global_conf->receiver_z(std::stod(wp[3]));
	}
}

void CLI::update_config (const std::string& param) const
{
	if (!param.compare("--help"))
		this->global_conf->help_opt(true);

	if (!param.compare("--version"))
		this->global_conf->version_opt(true);

	if (!param.compare("--safe"))
		this->global_conf->safe_mode(true);
}

bool CLI::is_float ( const std::string& literal )
{
	std::istringstream iss(literal);
	float f;
	iss >> std::noskipws >> f; // noskipws considers leading whitespace invalid
	// Check the entire string was consumed and if either failbit or badbit is set
	return iss.eof() && !iss.fail(); 
}

bool CLI::is_int ( const std::string& literal )
{
	bool is_empty = literal.empty();
	bool isnt_corret_start = ((!isdigit(literal[0])) && (literal[0] != '-') && (literal[0] != '+'));
	if( is_empty || isnt_corret_start ) 
		return false;

	char* p;
	strtol(literal.c_str(), &p, 10);

	return (*p == 0);
}

bool CLI::is_pair (const std::string& literal)
{
	std::vector<std::string> pair = CLI::split(literal, ',');
	if (pair.size() != 2) return false;
	return CLI::is_float(pair[0]) && CLI::is_float(pair[1]);
}

bool CLI::is_path (const std::string& literal)
{
	std::ifstream test(literal); 
	if (test) return true;
	else return false;
}

bool CLI::is_const_value (const std::string& literal)
{
	if (literal.empty()) return false;
	std::size_t splitter = literal.find(':');

	if (splitter != std::string::npos) return false;
	return CLI::is_float(literal);
}

bool CLI::is_range_value (const std::string& literal)
{
	if (literal.empty()) return false;
	std::vector<std::string> ranged = CLI::split(literal, ':');
	if (ranged.size() != 3) return false;

	bool is_float = CLI::is_float(ranged[0]) && CLI::is_float(ranged[1]) && CLI::is_float(ranged[2]);
	
	if (is_float) return std::stod(ranged[0]) < std::stod(ranged[2]);
	else return false;
}

bool CLI::is_onedim_model (const std::string& literal)
{
	std::vector<std::string> world_point = CLI::split(literal, ',');
	if (world_point.size() != 4) return false;

	std::size_t ranges = 0;
	std::size_t consts = 0;

	for (auto param : world_point) {
		if (CLI::is_range_value(param)) ranges++;
		if (CLI::is_const_value(param)) consts++;
	}

	return (ranges == 1) && (ranges+consts == 4);
}

bool CLI::is_twodim_model (const std::string& literal)
{
	std::vector<std::string> world_point = CLI::split(literal, ',');

	if (world_point.size() != 4) return false;
	std::size_t ranges = 0;
	std::size_t consts = 0;

	for (auto param : world_point) {
		if (CLI::is_range_value(param)) ranges++;
		if (CLI::is_const_value(param)) consts++;
	}

	return (ranges == 2) && (ranges+consts == 4);
}

bool CLI::is_event (const std::string& literal)
{
	std::vector<std::string> world_point = CLI::split(literal, ',');
	if (world_point.size() != 4) return false;
	
	bool is_floats = true;

	for (auto i : world_point) {
		bool condition = CLI::is_float(i) || CLI::is_range_value(i);
		if (!condition)
			is_floats = false;
	}

	if (is_floats) {
		double phi = std::stod(world_point[2]);
		if (phi >= 0 && phi <= 360) return true;
	}

	return false;
}

std::vector<std::string> CLI::split (const std::string& literal, char separator)
{
	std::stringstream test(literal);
	std::string segment;
	std::vector<std::string> seglist;

	while(std::getline(test, segment, separator))
	   seglist.push_back(segment);

	return seglist;
}

void CLI::set_mvc_model(PlotModel* model)
{
	this->plot_model = model;
}

void CLI::print_version () const
{
	std::cout << "Maxwell 0.9 private beta" << std::endl;
	std::cout << "Copyright (C) 2017 Rolan Akhmedov (Ukraine)" << std::endl;
	std::cout << "Contact e-mail: rolan.kharkiv@gmail.com" << std::endl;
}

void CLI::print_help () const
{
	// [--renoise] 

	std::cout << "usage: maxwell [--version] [--help] [--safe] [--conf <posix_path>]" << std::endl;
	std::cout << "               [--shape <type>] [--duration <float>] [<MODEL> <model_num>]" << std::endl;
	std::cout << "               [--magnitude <float>] [--radius <float>] [--kerr <float>]" << std::endl;
	std::cout << "               [--mur <float>] [--epsr <float>] [--noise <percent,percent>]" << std::endl;
	std::cout << "               [<SYSTEM> <float:float:float,float,float,float>]" << std::endl;

	std::cout << std::endl;
	std::cout << "  <MODEL>         --plot                --data                          " << std::endl;
	std::cout << "------------------------------------------------------------------------" << std::endl;
	std::cout << "        1         Ex(ct)                Ex(phi, rho, z)                 " << std::endl;
	std::cout << "        2         Hy(ct)                Hy(phi, rho, z)                 " << std::endl;
	std::cout << "        3         Ex(ct,rho)                                            " << std::endl;
	std::cout << "        4         Hy(ct,rho)                                            " << std::endl;
	std::cout << "        5         Ex(ct,z)                                              " << std::endl;
	std::cout << "        6         Hy(ct,z)                                              " << std::endl;
	std::cout << "        7         Hy(ct,N=1e1,1e2,1e3)                                  " << std::endl;

	std::cout << std::endl;
	std::cout << "  <SYSTEM>        --cartesian           --cylindric                     " << std::endl;
	std::cout << "------------------------------------------------------------------------" << std::endl;
	std::cout << "      f(ct)       <from>:<step>:<to>,<value>,<value>,<value>            " << std::endl;
	std::cout << "  f(ct,phi)       <from>:<step>:<to>,<float>,<from>:<step>:<to>,<float> " << std::endl;

	std::cout << std::endl;
	std::cout << "  --shape         DEFAULT_DURATION     OTHER_PARAMITERS                 " << std::endl;
	std::cout << "------------------------------------------------------------------------" << std::endl;
	std::cout << "      rect                infinity                                      " << std::endl;
	std::cout << "     gauss                       1              sigma=1                 " << std::endl;
	std::cout << "  triangle                       1                                      " << std::endl;

	std::cout << std::endl;
	std::cout << "  PARAMITER         DEFAULT_VALUE      OTHER_OPTIONS                    " << std::endl;
	std::cout << "------------------------------------------------------------------------" << std::endl;	
	std::cout << "  --magnitude                   1      any positive float               " << std::endl;
	std::cout << "  --radius                      1      any positive float               " << std::endl;
	std::cout << "  --mur                         1      any positive float               " << std::endl;
	std::cout << "  --epsr                        1      any positive float               " << std::endl;
	std::cout << "  --shape                    rect      gauss,triangle                   " << std::endl;
	std::cout << "  --noise                       0      any positive float (%)           " << std::endl;
	std::cout << "  --kerr                        0      any positive float               " << std::endl;


	std::cout << std::endl;
	std::cout << "All numerical options are equal to zero or grater" << std::endl;
	std::cout << "More options are avalible in maxwell.conf" << std::endl;
}

void CLI::print_arguments () const
{
	std::cout << std::boolalpha;
	std::cout << "Configuration of runtime..." << std::endl;
	
	auto t = this->global_conf->receiver_vt(); 
	auto x = this->global_conf->receiver_rho(); 
	auto y = this->global_conf->receiver_phi(); 
	auto z = this->global_conf->receiver_z();
	std::cout << "World point (m): ";
	if (t[0] == t[2]) std::cout << t[0] << ", ";
	else std::cout << t[0] << ':' << t[1] << ':' << t[2] << ", ";
	if (x[0] == x[2]) std::cout << x[0] << ", ";
	else std::cout << x[0] << ':' << x[1] << ':' << x[2] << ", ";
	if (y[0] == y[2]) std::cout << y[0] << ", ";
	else std::cout << y[0] << ':' << y[1] << ':' << y[2] << ", ";
	if (z[0] == z[2]) std::cout << z[0] << std::endl;
	else std::cout << z[0] << ':' << z[1] << ':' << z[2] << std::endl;

	std::cout << "Source radius (m): " << this->global_conf->plane_disk_radius() << std::endl;
	std::cout << "Source magnitude (m): " << this->global_conf->plane_disk_magnitude() << std::endl;
	std::cout << "Relative epsilon: " << this->global_conf->plane_disk_epsr() << std::endl;
	std::cout << "Relative mu: " << this->global_conf->plane_disk_mur() << std::endl;
	std::cout << "Kerr medium coefficient: " << this->global_conf->kerr_value() << std::endl;
	std::cout << "Conductivity (sigma): " << 0 << std::endl;
	
	#ifdef AUW_NOICE
		std::string noise_level = "Noise level: type(WHITE_UNIFORM), power(0), mu(FIRST), sigma(SECOND)";
	#else
		std::string noise_level = "Noise level: type(WHITE_GAUSS), power(0) mu(FIRST), sigma(SECOND)";
	#endif
	noise_level = noise_level.replace( noise_level.find("FIRST"), 5, std::to_string(0));
	std::string noise_sigma = std::to_string(this->global_conf->noise_level());
	noise_level = noise_level.replace( noise_level.find("SECOND"), 6, noise_sigma);
	std::cout << noise_level << std::endl;

	std::cout << "Magnetic component terms number: " << this->global_conf->magnetic_term_num() << std::endl;
	std::cout << "Float bitrae of GMP: " << this->global_conf->float_bitrate() << std::endl;
	std::cout << "Calculation thread number: " << this->global_conf->thread_number() << std::endl;
	
	if (this->global_conf->safe_mode()) {
		std::string server = this->global_conf->mysql_username() + "@" + this->global_conf->mysql_hostname(); 
		std::cout << "MySQL client will try to connect to " << server << "..." << std::endl;
	} else {
		std::cout << "No data base conection." << std::endl;
		std::cout << "Use --safe flag to store progress on hard drive." << std::endl;
	}
	std::cout << "..." << std::endl;
}

void CLI::read_config_file () const
{
	std::ifstream file(this->global_conf->maxwell_config_path());
	std::string line;

	while (std::getline(file, line)) {

		std::string new_line;

		unique_copy (line.begin(), line.end(), 
			std::back_insert_iterator<std::string>(new_line),
			[] (char a, char b) { return std::isspace(a) && std::isspace(b); 
		});

		if (new_line[0] == ' ') new_line.erase(0,1);
		std::size_t lsize = new_line.size() - 1;
		if (new_line[lsize] == ' ') new_line.erase(lsize);

		if (new_line.find("#") == std::string::npos) {
			std::size_t eq_char = new_line.find("=");
			if (eq_char != std::string::npos) {

				std::string arg = new_line.substr(eq_char+1);
				if (arg[0] == ' ') arg.erase(0,1);
				std::size_t arg_size = arg.size() - 1;
				if (arg[arg_size] == ' ') arg.erase(arg_size);
				std::string option = new_line.substr(0, eq_char);

				if (option.find("PRINT_PROGRESS") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						this->global_conf->print_progress(true);
					else if (arg.find("FALSE") != std::string::npos)
						this->global_conf->print_progress(false);
					else throw std::logic_error("Only TRUE/FALSE allowed...");
				}

				else if (option.find("DISPLAY_PARAMS") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						this->global_conf->display_params(true);
					else if (arg.find("FALSE") != std::string::npos)
						this->global_conf->display_params(false);
					else throw std::logic_error("Only TRUE/FALSE allowed...");
				}

				else if (option.find("CALL_GNUPLOT_BIN") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						this->global_conf->call_gnuplot(true);
					else if (arg.find("FALSE") != std::string::npos)
						this->global_conf->call_gnuplot(false);
					else throw std::logic_error("Only TRUE/FALSE allowed...");
				}

				else if (option.find("PLOT_COLOR_MAP") != std::string::npos) {
					if (arg.find("GRAYSCALE") != std::string::npos)
						this->global_conf->plot_color_map(Colormap::gray);
					else if (arg.find("PARULA") != std::string::npos)
						this->global_conf->plot_color_map(Colormap::parula);
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("PLOT_GRID") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						this->global_conf->plot_grid(true);
					else if (arg.find("FALSE") != std::string::npos)
						this->global_conf->plot_grid(false);
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("PLOT_BOUND_CAGE") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						this->global_conf->plot_baund_cage(true);
					else if (arg.find("FALSE") != std::string::npos)
						this->global_conf->plot_baund_cage(false);
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("PLOT_DEVICE") != std::string::npos) {
					if (arg.find("X11") != std::string::npos)
						this->global_conf->plot_device(PlotDev::x11);
					else if (arg.find("TERM") != std::string::npos)
						this->global_conf->plot_device(PlotDev::term);
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("GNUPLOT_BIN") != std::string::npos) {
					// std::experimental::filesystem::exists(arg);
					this->global_conf->path_gnuplot_binary(arg);
				}

				else if (option.find("SCRIPT_NAME") != std::string::npos) {
					// std::experimental::filesystem::exists(arg);
					this->global_conf->gnp_script_path(arg);
				}

				else if (option.find("MAGNETIC_TERM_NUM") != std::string::npos) {
					this->global_conf->magnetic_term_num(std::stoi(arg));
				}

				else if (option.find("FLOAT_BITRATE") != std::string::npos) {
					this->global_conf->float_bitrate(std::stoi(arg));
				}

				else if (option.find("PDISK_RADIUS") != std::string::npos) {
					this->global_conf->plane_disk_radius(std::stof(arg));
				}

				else if (option.find("PDISK_MAGNITUDE") != std::string::npos) {
					this->global_conf->plane_disk_magnitude(std::stof(arg));
				}

				else if (option.find("PDISK_EPSR") != std::string::npos) {
					this->global_conf->plane_disk_epsr(std::stof(arg));
				}

				else if (option.find("PDISK_MUR") != std::string::npos) {
					this->global_conf->plane_disk_mur(std::stof(arg));
				}

				else if (option.find("THREAD_NUMBER") != std::string::npos) {
					this->global_conf->thread_number(std::stof(arg));
				}

				else if (option.find("SUPERPOSITION") != std::string::npos) {
					if (arg.find("ADDITIVE") != std::string::npos) {
						this->global_conf->medium_superposition(Superposition::additive);
					} else if (arg.find("MULTIPL") != std::string::npos) {
						this->global_conf->medium_superposition(Superposition::multipl);
					} else {
						std::string text = "The argument is not allowed of SUPERPOSITION";
						throw std::logic_error(text);
					}
				}

				else if (option.find("MYSQL_HOSTNAME") != std::string::npos) {
					this->global_conf->mysql_hostname(arg);
				}

				else if (option.find("MYSQL_USERNAME") != std::string::npos) {
					this->global_conf->mysql_username(arg);
				}

				else if (option.find("MYSQL_PASSWORD") != std::string::npos) {
					this->global_conf->mysql_password(arg);
				}

				else if (option.find("MYSQL_DATABASE") != std::string::npos) {
					this->global_conf->mysql_database(arg);
				}
			}
		} 
	}
	file.close();
}
