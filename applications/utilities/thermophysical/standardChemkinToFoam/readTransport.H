std::vector<std::vector<std::string> > V;
std::ifstream ifs(args[3]);

std::string str;

while (getline(ifs,str))
{
	std::vector<std::string> vs;

	std::istringstream is(str);

	std::string s;

	while (is >> s)
	{
		if (s != "!")
		{
			vs.push_back(s);
		}
		else
		{
			break;
		}
	}

	V.push_back(vs);
}

std::unordered_map<std::string, std::vector<double> > m;

for (long unsigned int i = 0; i < V.size(); ++i)
{
	std::string name;
	std::vector<double> v;

	for (long unsigned int j = 0; j < V[i].size(); ++j)
	{
		if (j == 0)
		{
			name = V[i][j];
		}
		else
		{
			v.push_back(atof((V[i][j]).c_str()));
		}
	}
	m.insert(make_pair(name,v));
}

/*
for (std::unordered_map<std::string, std::vector<double> >::iterator it = m.begin(); it != m.end(); ++it)
{
	Info << it->first << " ";
	for (std::vector<double>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
	{
		Info << *jt << " ";
	}
	Info << endl;
}
*/


/*
dictionary specDict;

for (std::unordered_map<std::string, std::vector<double> >::iterator it = m.begin(); it != m.end(); ++it)
{
	dictionary dict("transport");

    dict.add("geometry", (it->second)[0]);
    dict.add("welldepth", (it->second)[1]);
    dict.add("diameter", (it->second)[2]);
    dict.add("dipole", (it->second)[3]);
    dict.add("polar", (it->second)[4]);
    dict.add("relaxation", (it->second)[5]);

	dictionary transDict;

	transDict.add("transport", dict);

	word name = it->first;
	specDict.add(name, transDict);
}

OFstream ofs("transport.dat");

specDict.write(ofs, false);
*/


OFstream ofs("transport.dat");

for (std::unordered_map<std::string, std::vector<double> >::iterator it = m.begin(); it != m.end(); ++it)
{
	word name = it->first;

	ofs.beginBlock(name);

	// Entries in dictionary format
	{
	    ofs.beginBlock("transport");
	    ofs.writeEntry("geometry", (it->second)[0]);
	    ofs.writeEntry("welldepth", (it->second)[1]);
	    ofs.writeEntry("diameter", (it->second)[2]);
	    ofs.writeEntry("dipole", (it->second)[3]);
	    ofs.writeEntry("polar", (it->second)[4]);
	    ofs.writeEntry("relaxation", (it->second)[5]);
	    ofs.endBlock();
	}

	ofs.endBlock();

	ofs << endl;
}
