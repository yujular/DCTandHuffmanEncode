/*
*参考资料:
*DCT矩阵运算:https://www.cnblogs.com/wyuzl/p/7880124.html
*图片大小2160*3840
*/
#include<iostream>
#include<cmath>
#include<map>
#include<string>
#include<cstring>
#include<opencv2/opencv.hpp>
#include"HuffmanEncode.h"

struct IHnode {
	IHnode *lchild;
	IHnode *rchild;
	IHnode *parent;
	int ch;
	IHnode() :lchild(NULL), rchild(NULL), parent(NULL) {}
};
const int N = 8;//分块大小
const double Q1[8][8] = {
{16,11,10,16,24,40,51,61},
{12,12,14,19,26,58,60,55},
{14,13,16,24,40,57,69,56},
{14,17,22,29,51,87,80,62},
{18,22,37,56,68,109,103,77},
{24,35,55,64,81,104,113,92},
{49,64,78,87,103,121,120,101},
{79,92,95,98,112,100,103,99}
};
const double Q2[8][8] = {
{17,18,24,47,99,99,99,99},
{18,21,26,66,99,99,99,99},
{24,26,56,99,99,99,99,99},
{47,66,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
{99,99,99,99,99,99,99,99},
};
const int dct_orderItem[8][8] = {
{0,1,5,6,14,15,27,28},
{2,4,7,13,16,26,29,42},
{3,8,12,17,25,30,41,43},
{9,11,18,24,31,40,44,53},
{10,19,23,32,39,45,52,54},
{20,22,33,38,46,51,55,60},
{21,34,37,47,50,56,59,61},
{35,36,48,49,57,58,62,63}
};
const double pi = 3.1415926535;

using namespace std;
using namespace cv;

Hnode dct_node[10000000];//压缩时节点序列
ChCode dct_code[10000000];//压缩时编码序列
IHnode *idct_Hroot = new IHnode();//霍夫曼解码树根节点

int ilz_idx = 0;
int ilz_number = 0;
int ilz_maxNum = 5000000;

int dct_Total = 0;
int idct_Total = 0;
int irows=1080;
int icols=1920;


double dct_c(int x) {
	return x == 0 ? sqrt(1.0 / double(N)) : sqrt(2.0 / (double)(N));
}

void dct_getCoe(Mat &m) {//计算DCT系数矩阵
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			m.ptr<double>(i)[j] = dct_c(i)*cos((j + 0.5)*pi*i / N);
		}
}

void dct_dealimg(const Mat &coe, Mat &img, const Mat &Q) {//对图片矩阵img进行DCT和量化
	for (int i = 0; i < img.rows; i += 8) {
		for (int j = 0; j < img.cols; j += 8) {
			img(Range(i, i + 8), Range(j, j + 8)) = coe * img(Range(i, i + 8), Range(j, j + 8))*coe.t();//DCT,F=AfAT
			img(Range(i, i + 8), Range(j, j + 8)) /= Q;
		}
	}
}

void dct_idealimg(const Mat &coe, Mat&img, const Mat &Q) {//对图片反量化，反DCT
	for (int i = 0; i < img.rows; i += 8) {
		for (int j = 0; j < img.cols; j += 8) {
			img(Range(i, i + 8), Range(j, j + 8)) = img(Range(i, i + 8), Range(j, j + 8)).mul(Q);//反量化
			img(Range(i, i + 8), Range(j, j + 8)) = coe.t()*img(Range(i, i + 8), Range(j, j + 8))*coe;//反DCT
		}
	}
}

void getHuffmanSeq(Mat &img) {
	ofstream fout("HuffmanCodeSeq.bin", ios::out);//打开二进制文件
	ofstream fout2("HuffmanCodeSeq.txt",ios::out);//记录01序列
	ofstream fout3("SeqLength.txt", ios::out || ios::app);//记录序列长度
	long long dct_length = 0;
	unsigned char buf = 0;//位输出
	int count = 0;//8位计数器
	std::map <int, deque<int>> m;//码表
	for (int i = 0; i < dct_Total; i++) {
		m[dct_code[i].ch] = dct_code[i].hoffmanCode;
	}
	for (int i = 0; i < img.rows; i += 8) {//遍历图片
		for (int j = 0; j < img.cols; j += 8) {
			deque<int> dct_seq[64];
			for (int k = 0; k < 8; k++) {//遍历图片
				for (int l = 0; l < 8; l++) {
					int now_pos = dct_orderItem[k][l];//当前元素位置
					int now_char = img.ptr<short>(i + k)[j + l];
					dct_seq[now_pos] = m[now_char];
				}
			}
			for (int m = 0; m < 64; m++) {//输出当前64个序列
				deque<int>::iterator iter;
				dct_length += dct_seq[m].size();
				for (iter = dct_seq[m].begin(); iter != dct_seq[m].end(); iter++) {
					buf = buf | (*iter) << (7 - count);
					fout2 << *iter;
					count++;
					if (count >= 8) {
						count = 0;
						fout << buf;
						buf = 0;
					}
				}
			}
		}
	}
	if (count != 0) {
		fout << buf;
	}
	fout3 << dct_length << endl;
	fout.close();
	fout2.close();
	fout3.close();
}

void getHuffmanSeq2(Mat &img) {
	ofstream fout3("SeqLength2.txt", ios::out || ios::app);//记录序列长度
	long long dct_length = 0;
	std::map <int, deque<int>> m;//码表
	for (int i = 0; i < dct_Total; i++) {
		m[dct_code[i].ch] = dct_code[i].hoffmanCode;
	}
	for (int i = 0; i < img.rows; i += 8) {//遍历图片
		for (int j = 0; j < img.cols; j += 8) {
			deque<int> dct_seq[64];
			for (int k = 0; k < 8; k++) {//遍历图片
				for (int l = 0; l < 8; l++) {
					int now_pos = dct_orderItem[k][l];//当前元素位置
					int now_char = img.ptr<short>(i + k)[j + l];
					dct_seq[now_pos] = m[now_char];
				}
			}
			for (int m = 0; m < 64; m++) {//输出当前64个序列
				dct_length += dct_seq[m].size();
			}
		}
	}
	fout3 << dct_length << endl;
	fout3.close();
}

void RebulidHuffmanTree()
{
	ifstream fin("HuffmanCodeTable.txt", ios::in);
	int n;
	fin >> n;
	idct_Total = n;
	for (int i = 0; i < n; i++)
	{
		int now_ch;
		fin >> now_ch;//当前符号
		IHnode *now_node = idct_Hroot;//当前节点
		char buf[256] = {0};//最长编码长度256以内
		fin.getline(buf, 256);
		for (int j = 0; buf[j] != 0; j++) {
			if (buf[j] == '0') {
				if (now_node->lchild == NULL) {
					now_node->lchild = new IHnode();//创建新节点
					now_node->lchild->parent = now_node;//父节点赋值
				}
				now_node = now_node->lchild;
			}
			else if (buf[j] == '1') {
				if (now_node->rchild == NULL) {
					now_node->rchild = new IHnode();//创建新节点
					now_node->rchild->parent = now_node;//父节点赋值
				}
				now_node = now_node->rchild;
			}
		}
		now_node->ch = now_ch;
	}
	fin.close();
}

void dct_huffmanEncode(Mat &img) {
	map<int, int> Map;
	ofstream fout("pro1.txt", ios::out);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			int now = img.ptr<short>(i)[j];
			Map[now]++;//统计出现次数
		}
	}
	map<int, int>::iterator iter;
	for (iter = Map.begin(); iter != Map.end(); iter++) {
		dct_code[dct_Total].ch = iter->first;
		dct_node[dct_Total].weight = iter->second;
		fout << iter->first << "  " << iter->second << endl;
		dct_Total++;
	}
	priority_queue<Hnode*, vector<Hnode*>, cmp> pq;//优先队列
	for (int i = 0; i < dct_Total; i++)//节点入队列
		pq.push(dct_node + i);
	getHuffmanCode(pq, dct_code, dct_node, dct_Total);//求霍夫曼编码
	getHuffmanTable(dct_code, dct_Total);//输出霍夫曼码表
	getHuffmanSeq(img);//输出霍夫曼序列
}
void getHuffmanTable2(ChCode *code, int n)
{
	ofstream fout("HuffmanCodeTable2.txt", ios::out);
	fout << n << endl;
	for (int i = 0; i < n; i++)
	{
		fout << code[i].ch << " ";
		deque<int>::iterator iter;
		for (iter = code[i].hoffmanCode.begin(); iter != code[i].hoffmanCode.end(); iter++)
			fout << *iter;
		fout << endl;
	}
	fout.close();
}
void dct_huffmanEncode2(Mat &img) {
	map<int, int> Map;
	ofstream fout("pro2.txt", ios::out);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			int now = img.ptr<short>(i)[j];
			Map[now]++;//统计出现次数
		}
	}
	map<int, int>::iterator iter;
	for (iter = Map.begin(); iter != Map.end(); iter++) {
		dct_code[dct_Total].ch = iter->first;
		dct_node[dct_Total].weight = iter->second;
		fout << iter->first << "  " << iter->second << endl;
		dct_Total++;
	}
	priority_queue<Hnode*, vector<Hnode*>, cmp> pq;//优先队列
	for (int i = 0; i < dct_Total; i++)//节点入队列
		pq.push(dct_node + i);
	getHuffmanCode(pq, dct_code, dct_node, dct_Total);//求霍夫曼编码
	getHuffmanTable2(dct_code, dct_Total);//输出霍夫曼码表
	getHuffmanSeq2(img);//输出霍夫曼序列
}


void dct_huffmanDecode(Mat &img) {
	RebulidHuffmanTree();//重建霍夫曼树
	ifstream fin("HuffmanCodeSeq.bin", ios::in);//读入序列 
	deque<int> deq;
	for (int i = 0; i < irows; i += 8) {
		for (int j = 0; j < icols; j += 8) {
			int ch[64];//记录当前矩阵下每个元素的值
			for (int k = 0; k < 64; k++) {
				IHnode *now_node = idct_Hroot;//从根节点寻找
				while (1) {
					if (deq.empty()) {
						unsigned char cch;
						cch = fin.get();
						if (cch != '\r') {
							stringstream a;
							for (int x = 7; x > 0; x--) {
								a << ((cch >> x) & 1);
							}
							a << (cch & 1);
							string s = a.str();
							for (int x = 0; x < 8; x++) {
								deq.push_back(s[x] - '0');
							}
						}
					}
					int now_code = deq.front();//取第一个码字
					deq.pop_front();//出队列
					if (now_code == 1) {
						if (now_node->rchild != NULL) {//可继续dfs
							now_node = now_node->rchild;
						}
						else {//不可继续dfs
							deq.push_front(now_code);//当前码字属于下一个，归还队列
							ch[k] = now_node->ch;//找到当前字符
							now_node = idct_Hroot;
							break;
						}
					}
					else if (now_code == 0) {
						if (now_node->lchild != NULL) {//可继续dfs
							now_node = now_node->lchild;
						}
						else {//不可继续dfs
							deq.push_front(now_code);//当前码字属于下一个，归还队列
							ch[k] = now_node->ch;//找到当前字符
							now_node = idct_Hroot;
							break;
						}
					}
				}
			}
			for (int k = 0; k < 8; k++) {
				for (int l = 0; l < 8; l++) {
					int t = ch[dct_orderItem[k][l]];
					img.ptr<short>(i + k)[j + l] = t;
				}
			}
		}
	}
}

void dct_getQ(int q[8][8],int x) {
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			q[i][j] = 1 + x * (i + j - 1);
		}
	}
}

class TreeNode {
public:
	int index;
	char entry;
	TreeNode *parent;
	TreeNode *firstChild;
	TreeNode *sibling;

	TreeNode()
	{
		entry = NULL;
		index = NULL;
		parent = NULL;
		firstChild = NULL;
		sibling = NULL;
	};
	TreeNode(char e, int inx, TreeNode *par)
	{
		entry = e;
		index = inx;
		par = parent;
		firstChild = NULL;
		sibling = NULL;
	};
	TreeNode* findChild(char c)
	{
		if (firstChild == NULL) return NULL;
		TreeNode *kid = firstChild;
		while (kid != NULL)
		{
			if (kid->entry == c) return kid;
			kid = kid->sibling;
		}
		return NULL;
	};
	void addChild(char c, int mycount)
	{
		if (firstChild == NULL)
		{
			firstChild = new TreeNode(c, mycount, this);
		}
		else
		{
			TreeNode *kid = firstChild;
			while (kid->sibling != NULL) kid = kid->sibling;
			kid->sibling = new TreeNode(c, mycount, this);
		}
	};

};

void LZencode(string s)
{
	ofstream outfile;
	outfile.open("LZ78encoding.txt");
	int total = s.length();
	char c;
	TreeNode *root = new TreeNode(NULL, 0, NULL);
	TreeNode *cur = root;
	TreeNode *loc = cur;
	for (int i = 0; i < total; i++)
	{
		c = s[i];
		loc = cur->findChild(c);

		if (loc == NULL)
		{
			cur->addChild(c, i);
			outfile << cur->index << c << " ";
			cur = root;
		}
		else { cur = loc; }
	}
	outfile.close();
}


void dct_lzencode(Mat &img) {
	string s;
	for (int i = 0; i < img.rows; i += 8) {
		for (int j = 0; j < img.cols; j += 8) {
			char ch[64];
			for (int k = 0; k < 8; k++) {
				for (int l = 0; l < 8; l++) {
					ch[dct_orderItem[k][l]] = img.ptr<short>(i+k)[j+l]+100;
				}
			}
			string ss = ch;
			s += ss;
		}
	}
	LZencode(s);//lz编码
}
struct ilz_Node {
	int index;
	char *code;
};
ilz_Node *ilz_List;
void doubleSpace()
{
	ilz_Node *tmp = ilz_List;
	ilz_List = new ilz_Node[2 * ilz_maxNum];
	for (int i = 0; i < ilz_maxNum; i++)
	{
		ilz_List[i].index = tmp[i].index;
		ilz_List[i].code = new char[strlen(tmp[i].code) + 1];
		strcpy(ilz_List[i].code, tmp[i].code);
	}
}
int str_to_num(char *str, int len)
{
	int ans = 0;
	int base = 1;
	while (len)
	{
		ans += (str[len - 1] - '0')*base;
		len--;
		base *= 10;
	}
	return ans;
}
void lzdecode() {
	ilz_List = new ilz_Node[ilz_maxNum];
	ifstream infile("LZ78Encoding.txt");
	ofstream outfile("LZ78decoding.txt");
	char ch;
	if (infile.fail()) cout << "open error" << endl;
	infile.get(ch);
	while (ch != '@')
	{
		ilz_number++;
		//if(number>maxNum) doubleSpace();
		char tmpstr[20];
		int tmplen = 0;
		while (ch != ' ' && ch != '@')
		{
			if (ch != '\r')//过滤掉Win系统会产生的回车
			{
				tmpstr[tmplen] = ch;
				tmplen++;
			}
			if (ch != '@')infile.get(ch);
		}
		//空格有两种情况，一种是编码中的空格，一种是分隔符，这时我们要根据空格后面是不是还是空格来判断是不是分隔符
		if (ch == ' ') infile.get(ch);//读空格后面的一个字符
		if (ch == ' ')//空格后面还是空格，证明上一个空格应该被写进去，但是这个空格是分隔符所以不要了
		{
			tmpstr[tmplen] = ' ';
			tmplen++;
			infile.get(ch);
		}//剩下一种情况分别有数字和EOF两种情况，我们就交给下一轮处理

		ilz_List[ilz_idx].index = str_to_num(tmpstr, tmplen - 1) - 1;//如果是0设成-1以此类推
		if (ilz_List[ilz_idx].index == -1)
		{
			ilz_List[ilz_idx].code = new char[2];
			ilz_List[ilz_idx].code[0] = tmpstr[tmplen - 1];
			ilz_List[ilz_idx].code[1] = '\0';
		}
		else {
			int codelen = strlen(ilz_List[ilz_List[ilz_idx].index].code);
			ilz_List[ilz_idx].code = new char[codelen + 2];
			strcpy(ilz_List[ilz_idx].code, ilz_List[ilz_List[ilz_idx].index].code);
			ilz_List[ilz_idx].code[codelen] = tmpstr[tmplen - 1];
			ilz_List[ilz_idx].code[codelen + 1] = '\0';
		}
		outfile << (ilz_List[ilz_idx].code)<<" ";
		ilz_idx++;
	}
	infile.close();
	outfile.close();
}
int main() {
	Mat dct_coefficient(Size(8, 8), CV_64FC1);//DCT系数矩阵
	Mat dct_Q1(Size(8, 8), CV_64FC1, (Scalar*)Q1);//量化系数矩阵
	Mat dct_Q2(Size(8, 8), CV_64FC1, (Scalar*)Q2);//量化系数矩阵2
	/*******************************************分割线***********************************************/
	/*******************************************压缩***********************************************/
	Mat img;
	img = imread("img.jpg", CV_LOAD_IMAGE_GRAYSCALE);//插入图片,灰度图
	if (img.empty()) {
		cout << "打开图片失败,请检查路径" << endl;
		return 0;
	}
	imwrite("beforedeal.jpg", img);
	Mat dealimg;//转换成double矩阵
	img.convertTo(dealimg, CV_64FC1);
	dct_getCoe(dct_coefficient);//获取DCT系数矩阵
	dct_dealimg(dct_coefficient, dealimg, dct_Q1);//DCT变化，量化
	Mat intimg;
	dealimg.convertTo(intimg, CV_16SC1);//对量化后图片取整
	//dct_lzencode(intimg);//lz78编码
	//dct_huffmanEncode(intimg);//霍夫曼编码
	/*******************************************分割线***********************************************/
	/*******************************************解压缩***********************************************/
	Mat i_intimg(Size(icols,irows),CV_16SC1);
	Mat i_dealimg;
	Mat i_img;
	//dct_huffmanDecode(i_intimg);//从霍夫曼编码中获取图片矩阵
	intimg.convertTo(i_dealimg, CV_64FC1);//转成double矩阵
	dct_idealimg(dct_coefficient, i_dealimg, dct_Q1);//DCT反量化，反变化
	i_dealimg.convertTo(i_img, CV_8UC1);//取整输出
	cv::imwrite("afterdeal.jpg", i_img);//输出图片
	/*******************************************分割线***********************************************/
	/****************************************不同量化矩阵*********************************************/
	Mat dealimg2;
	img.convertTo(dealimg2, CV_64FC1);
	cout << dealimg2(Range(0, 16), Range(0, 16)) << endl;
	dct_dealimg(dct_coefficient, dealimg2, dct_Q2);
	Mat initimg2;
	Mat i_img2;
	dealimg2.convertTo(initimg2, CV_16SC1);
	cout << initimg2(Range(0, 16), Range(0, 16)) << endl;
	//dct_huffmanEncode2(initimg2);
	dct_lzencode(initimg2);
	initimg2.convertTo(dealimg2, CV_64FC1);//转成double矩阵
	dct_idealimg(dct_coefficient, dealimg2, dct_Q2);//DCT反量化，反变化
	dealimg2.convertTo(i_img2, CV_8UC1);//取整输出
	cv::imwrite("afterdeal2.jpg", i_img2);//输出图片

	return 0;
}
