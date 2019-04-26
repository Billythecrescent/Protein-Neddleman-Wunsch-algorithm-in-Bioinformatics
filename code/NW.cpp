//头文件 
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;
#define BLANK -8

//定义NW表单位的结构体
struct Unit {
	int value;
	string p; //p表示回溯时的路线。为diag表示走对角线，p为left表示走左，p为up表示向上
};

//读入蛋白质序列
//string protein1 = "IPGAWD";
//string protein2 = "VGAWAD";


int max(int a, int b, int c) {
	return (a >= b ? a : b) >= c ? (a >= b ? a : b) : c;
}

//判断一个char是否在vector中
bool is_element_in_vector(vector<char> v, char element) {
	vector<char>::iterator it;
	it = find(v.begin(), v.end(), element);
	if (it != v.end()) {
		return true;
	}
	else{
		return false;
	}
}

//存储BLOSUM氨基酸对照表

/*void storeBLOSUMtable(map<string, int>& BLOSUM) {

	BLOSUM["VI"] = 3;
	BLOSUM["VP"] = -3;
	BLOSUM["VG"] = -3;
	BLOSUM["VA"] = 0;
	BLOSUM["VW"] = -3;
	BLOSUM["VD"] = -3;
	BLOSUM["GI"] = -4;
	BLOSUM["GP"] = -2;
	BLOSUM["GG"] = 7;
	BLOSUM["GA"] = 0;
	BLOSUM["GW"] = -2;
	BLOSUM["GD"] = -1;
	BLOSUM["WI"] = -2;
	BLOSUM["WP"] = -3;
	BLOSUM["WG"] = -2;
	BLOSUM["WA"] = -2;
	BLOSUM["WW"] = 15;
	BLOSUM["WD"] = -4;
	BLOSUM["AI"] = -1;
	BLOSUM["AP"] = -1;
	BLOSUM["AG"] = 0;
	BLOSUM["AA"] = 5;
	BLOSUM["AW"] = -2;
	BLOSUM["AD"] = -2;
	BLOSUM["DI"] = -4;
	BLOSUM["DP"] = -1;
	BLOSUM["DG"] = -1;
	BLOSUM["DA"] = -2;
	BLOSUM["DW"] = -4;
	BLOSUM["DD"] = 7;

}*/
vector<char> storeBLOSUMtable(const char* file, map<string, int>& BLOSUM) {
	string buf;
	ifstream fd(file, ios::in);
	vector<char> symbol;
	if (!fd)
	{
		cout << "Error in Opening " << file << endl;
		return symbol;
	}
	int n = 0;
	while (fd.good())
	{
		buf.clear();
		getline(fd, buf);

		//读取BLOSUM矩阵信息
		if (buf[0] == '#') {
			continue;
		}
		if (buf.length() > 0 && buf[0] == ' ') {
			istringstream stream(buf);
			while (stream) {
				char temp = ' ';
				stream >> temp;
				if (temp != ' ' && temp != '\n') {
					symbol.push_back(temp);
					n += 1;
				}
			}
		}
		if (buf.length() > 0 && buf[0] != ' ' && n != 0) {
			istringstream stream(buf);
			char str;
			stream >> str;
			//cout << str << endl;
			for (int i = 0; i < n; i++)
			{
				int value = INT_MIN;
				stream >> value;
				if (value != INT_MIN) {
					//cout << value << "  ";
					string index;
					index.push_back(str);
					index.push_back(symbol[i]);
					//cout << index << " ";
					BLOSUM[index] = value;
				}
			}
			//cout << endl;
		}
	}
	fd.close();

	return symbol;
}

//read the protein FASTA file
vector<string> readFASTA(const char* file) {
	string buf;
	vector<string> result;
	ifstream fd(file, ios::in);
	if (!fd)
	{
		cout << "Can not find the file " << file << endl;
		return result;
	}
	int n = 0;
	while (fd.good())
	{
		buf.clear();
		getline(fd, buf);
		if (buf[0] == '>') {
			string sequence;
			getline(fd, sequence);
			result.push_back(sequence);
		}
		else {
			if (buf != "" && buf != "\n" && buf != " ") {
				string temp = result.back();
				result.pop_back();
				temp += buf;
				result.push_back(temp);
			}
		}
	}
	return result;
}

//BLOSUM index
string BLOSUMindex(int i, int j, string sequence1, string sequence2) {
	string index;
	index.push_back(sequence1.at(i));
	index.push_back(sequence2.at(j));
	return index;
}

//构造NW算法表
//函数返回NW table
Unit** NWtable(string sequence1, string sequence2, map<string, int> BLOSUM) {
	const int n1 = sequence1.size(); 
	const int n2 = sequence2.size();
	Unit** table = new Unit*[n1+1]; // 分配所有行的首地址
	for (int i = 0; i < n1+1; i++) { // 按行分配每一列
		table[i] = new Unit[n2+1];
	}
	//table初始化 
	for (int i = 0; i < n1 + 1; i++) {
		for (int j = 0; j < n2 + 1; j++) {
			if (i == 0 && j == 0) {
				table[i][j].value = 0;
				table[i][j].p = "zero";
			}
			else if (i == 0 && j != 0) {
				table[i][j].value = j * BLANK;
				table[i][j].p = "left";
			}
			else if (i != 0 && j == 0) {
				table[i][j].value = i * BLANK;
				table[i][j].p = "upwa";
			}
			else {
				table[i][j].value = 0;
				table[i][j].p = "";
			}
		}
	}
	
	//填充table 
	for (int i = 1; i < n1+1; i++) {
		for (int j = 1; j < n2+1; j++) {
			string index = BLOSUMindex(i-1,j-1,sequence1,sequence2);
			/*cout << "BLOSUM" << endl;
			cout << i << " " << j << index << " " << BLOSUM[index] << endl;*/
			int maxvalue = max(table[i - 1][j - 1].value + BLOSUM[index], table[i-1][j].value+BLANK, table[i][j-1].value+BLANK);
			table[i][j].value = maxvalue;
			if (maxvalue == table[i - 1][j - 1].value + BLOSUM[index]) {
				table[i][j].p = "diag";
			}
			else if (maxvalue == table[i - 1][j].value + BLANK) {
				table[i][j].p = "upwa";
			}
			else if (maxvalue == table[i][j - 1].value + BLANK) {
				table[i][j].p = "left";
			}
			else {
				cout << "Error: Exceptional condition happens!" << endl;
			}
		}
	}

	return table;
}

//Unitback
int Unitback(vector<int> index, Unit**table, vector<vector<int> > &result) {
	result.push_back(index);
	if (table[index[0]][index[1]].p == "zero") {
		return 0;
	}
	else if (table[index[0]][index[1]].p == "diag") {
		vector<int> temp;
		temp.push_back(index[0] - 1);
		temp.push_back(index[1] - 1);
		return 1 + Unitback(temp, table, result);
	}
	else if (table[index[0]][index[1]].p == "left") {
		vector<int> temp;
		temp.push_back(index[0]);
		temp.push_back(index[1] - 1);
		return 0 + Unitback(temp, table, result);
	}
	else if (table[index[0]][index[1]].p == "upwa") {
		vector<int> temp;
		temp.push_back(index[0] - 1);
		temp.push_back(index[1]);
		return 0 + Unitback(temp, table, result);
	}
	else {
		return 0;
	}
}


//进行NW表的回溯计算
vector<vector<int> > NWtraceback(string sequence1, string sequence2,Unit** table) {
	const int n1 = sequence1.size();
	const int n2 = sequence2.size();
	vector<vector<int> > result;
	vector<int> index(2,-1);
	int maxvalue = INT_MIN;
	for (int i = 0; i < n1 + 1; i++) {
		if (table[i][n2].value > maxvalue) {
			maxvalue = table[i][n2].value;
			index[0] = i; index[1] = n2;
		}
	}
	for (int j = 0; j < n2 + 1; j++) {
		if (table[n1][j].value > maxvalue) {
			maxvalue = table[n1][j].value;
			index[0] = n1; index[1] = j;
		}
	}
	cout << "The max score is " << maxvalue << endl;
	int alignnum = Unitback(index, table, result);

	//result翻转就是结果
	reverse(result.begin(),result.end());

	//释放table空间
	for (int i = 0; i < n1 + 1; i++)
	{
		delete [] table[i];
	}
	delete [] table;
	cout << "Alignment proportion is " << alignnum*100/(n1 > n2 ? n2 : n1) << "%" << endl;
	return result;
}

//根据路径，输出匹配结果 
void Alignout(vector<vector<int> >result,string sequence1,string sequence2) {
	vector<char> seq1;
	vector<char> seq2;
	int n = result.size();
	for (int i = 1; i < n; i++) {
		if (result[i][0] > result[i - 1][0] && result[i][1] > result[i - 1][1]) {
			seq1.push_back(sequence1[result[i][0] - 1]);
			seq2.push_back(sequence2[result[i][1] - 1]);
		}
		else if (result[i][0] > result[i - 1][0] && result[i][1] == result[i - 1][1]) {
			seq1.push_back(sequence1[result[i][0] - 1]);
			seq2.push_back('-');
		}
		else if (result[i][0] == result[i - 1][0] && result[i][1] > result[i - 1][1]) {
			seq1.push_back('-');
			seq2.push_back(sequence2[result[i][1] - 1]);
		}
		else {
			cout << "Error: Exceptional condition happens!" << endl;
		}
	}
	
	//输出seq1，seq2
	for (int i = 0; i < n-1; i++) {
		cout << seq1[i];
	}
	cout << endl;
	for (int i = 0; i < n - 1; i++) {
		cout << seq2[i];
	}
	cout << endl;
}

int NWhandler(string sequence1, string sequence2, map<string, int> BLOSUM, vector<char>symbol) {
	int n1 = sequence1.size();
	int n2 = sequence1.size();
	
	cout << "Alignment between" << endl;
	cout << "sequence 1: " << sequence1 << endl;
	cout << "sequence 2: " << sequence2 << endl;
	
	/*char a = '*';
	cout << "a is in symbol? " << is_element_in_vector(symbol, a) << endl;*/
	for (int i = 0; i < sequence1.length(); i++) {
		if (!is_element_in_vector(symbol, sequence1[i])) {
			cout << "illegal character used in the first sequence." << endl;
			return 0;
		}
	}
	for (int i = 0; i < sequence2.length(); i++) {
		if (!is_element_in_vector(symbol, sequence2[i])) {
			cout << "illegal character used in the second sequence." << endl;
			return 0;
		}
	}
	cout << "Result: " << endl << "-----------------------------" << endl;;
	
	Unit** table = NWtable(sequence1, sequence2, BLOSUM);
	//测试table生成
	/*cout << "NW table is:" << endl;
	for (int i = 0; i < n1 + 1; i++)
	{
		for (int j = 0; j < n2 + 1; j++)
		{
			cout << table[i][j].value << " ";
		}
		cout << endl;
	}*/
	vector<vector<int> > trace = NWtraceback(sequence1, sequence2, table);
	Alignout(trace, sequence1, sequence2);
	return 1;
}


int main(int argc, char* argv[]) {
	//Needleman Wunsch Algoritm manupulation

	
	//string sequence1 = "VGAWADAGD"; 
	//string sequence2 = "IPGAWDWWI";

	if (*argv[1] == '-' &&  *(argv[1] + 1) == 'e') {
		string sequence1 = argv[2];
		string sequence2 = argv[3];
		vector<char> symbol;
		//得到BLOSUM45矩阵 
		map<string, int> BLOSUM;
		string filename = argv[4];
		string file = "../lib/" + filename;
		const char* path = file.data();

		symbol = storeBLOSUMtable(path, BLOSUM);
		cout << "Alignment using: " << filename << endl;
		if (symbol.size() == 0) {
			return 0;
		}

		if (NWhandler(sequence1, sequence2, BLOSUM, symbol)) {
			cout << "Alignment succeed" << endl << endl;
		}
		else {
			cout << "Alignment failed" << endl << endl;
			return 0;
		}
	}
	else if (*argv[1] == '-' &&  *(argv[1] + 1) == 'i') {
		string FASTA1 = argv[2];
		string file1 = "../sequence/" + FASTA1;
		const char* path1 = file1.data();
		string FASTA2 = argv[3];
		string file2 = "../sequence/" + FASTA2;
		const char* path2 = file2.data(); 
		vector<string> inputseq1 = readFASTA(path1);
		vector<string> inputseq2 = readFASTA(path2);
		/*for (int i = 0; i < inputseq2.size(); i++) {
			cout << inputseq2[i] << endl;
		}*/
		vector<char> symbol;
		//得到BLOSUM45矩阵 
		map<string, int> BLOSUM;
		string filename = argv[4];
		string file = "../lib/" + filename;
		const char* path = file.data();

		symbol = storeBLOSUMtable(path, BLOSUM);
		cout << "Alignment using: " << filename << endl;
		if (symbol.size() == 0) {
			return 0;
		}
		for (int i = 0; i < inputseq1.size(); i++) {
			for (int j = 0; j < inputseq2.size(); j++) {
				if (NWhandler(inputseq1[i], inputseq2[j], BLOSUM, symbol)) {
					cout << "Alignment succeed" << endl << endl;
				}
				else {
					cout << "Alignment failed" << endl << endl;
				}
			}
		}
	}
	return 0;
}
