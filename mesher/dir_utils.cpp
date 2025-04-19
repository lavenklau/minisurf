#define FMT_HEADER_ONLY
#include "fmt/core.h"
#include "dir_utils.h"
#include "regex"
#include <iostream>



std::vector<std::string> dir_utils::filterFiles(const std::vector<std::string>& files, std::string extname)
{
	std::vector<std::string> filteredfiles;
	for (int i = 0; i < files.size(); i++) {
		std::filesystem::path pth(files[i]);
		if (pth.extension().string() == extname) {
			filteredfiles.emplace_back(pth.string());
		}
	}
	return filteredfiles;
}

std::vector<std::string> dir_utils::matchFiles(const std::vector<std::string>& files, std::string regPattern)
{
	std::vector<std::string> matchedfiles;
	std::regex ptn(regPattern);
	for (int i = 0; i < files.size(); i++) {
		std::filesystem::path pth(files[i]);
		std::string fn = pth.filename().string();
		if (std::regex_match(fn, ptn)) {
			matchedfiles.emplace_back(files[i]);
		}
	}
	return matchedfiles;
}

std::vector<std::string> dir_utils::listFile(std::string dirname)
{
	std::filesystem::path pth(dirname);
	//std::filesystem::directory_iterator contents(pth);
	std::filesystem::recursive_directory_iterator contents(pth);

	std::vector<std::string> files;

	for (auto it : contents) {
		if (it.is_directory()) continue;
		
		files.emplace_back(it.path().string());
	}
	
	return files;
}

std::string dir_utils::path2filename(std::string pathstr)
{
	std::filesystem::path p(pathstr);
	return p.filename().string();
}

std::string dir_utils::path2extension(std::string pathstr)
{
	std::filesystem::path p(pathstr);
	return p.extension().string();
}

std::string dir_utils::incremental(std::string pathstr)
{
	std::filesystem::path p(pathstr);
	// if the path not exist, return the raw
	std::filesystem::path pat = p.parent_path();
	std::filesystem::path newPath = p;
	int id = 0;
	int iter = 0;
	while (std::filesystem::exists(newPath) || id == 0) {
		p = newPath;
		std::string fileName = p.stem().string();
		std::regex r("\\d+");
		std::smatch sma;
		bool flag = std::regex_search(fileName, sma, r);
		// fmt::print("sma[0] = {}\n", sma[0].str());
		if (flag) {
			// fmt::print("sma size = {}, stoi str = {}\n", sma.size(), sma[sma.size() - 1].str());
			id = std::stoi(sma[sma.size() - 1].str());
			std::string newid = fmt::format("{:04d}", ++id);
			fileName = std::regex_replace(fileName, r, newid.c_str());
			// fmt::print("new filename = {}\n", fileName);
		}
		else {
			std::string newid = fmt::format("{:04d}", id++);
			fileName = fileName + newid;
		}
		newPath = pat / (fileName + p.extension().string());
	}

	fmt::print("incremental file {}\n", newPath.string());
	return newPath.string();
	// exit(0);
}

bool dir_utils::exist(std::string path)
{
	std::filesystem::path p(path);
	return std::filesystem::exists(p);
}
