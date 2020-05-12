#include "HuffmanEncode.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <queue>

using namespace std;


void getHuffmanTree(priority_queue<Hnode *, vector<Hnode *>, cmp> &pq, int n)
{ //节点优先队列和节点名
	for (int i = 0; i < n; i++)
	{ //n次编码
		if (pq.size() == 1)
		{
			break; //只剩下一个节点
		}
		Hnode *new_Pnode = new Hnode();                            //新节点
		Hnode *new_Lnode = pq.top();                               //当前概率最小的
		pq.pop();                                                  //出队列
		Hnode *new_Rnode = pq.top();                               //第二小的
		pq.pop();                                                  //出队列
		new_Pnode->weight = new_Lnode->weight + new_Rnode->weight; //节点为两个的和
		new_Pnode->lchild = new_Lnode;
		new_Pnode->rchild = new_Rnode;
		new_Lnode->parent = new_Rnode->parent = new_Pnode;
		pq.push(new_Pnode);
	}
}

void getHuffmanCode(priority_queue<Hnode *, vector<Hnode *>, cmp> &pq, ChCode *code, Hnode *node, int n)
{
	getHuffmanTree(pq, n); //获取哈夫曼树
	for (int i = 0; i < n; i++)
	{
		Hnode *now_node = &node[i];
		while (now_node->parent != NULL)
		{
			if (now_node == now_node->parent->lchild)
			{
				code[i].hoffmanCode.push_front(0);
			}
			else
			{
				code[i].hoffmanCode.push_front(1);
			}
			now_node = now_node->parent;
		}
	}
}

void getHuffmanTable(ChCode *code, int n)
{
	ofstream fout("HuffmanCodeTable.txt", ios::out);
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
