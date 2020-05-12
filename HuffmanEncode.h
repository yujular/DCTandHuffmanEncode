
#pragma once

#ifndef _HuffmanEncode_H
#define _HuffmanEncode_H

#include<algorithm>
#include<queue>
#include<deque>
#include<iostream>
#include<vector>

typedef struct ChCode
{							//字符结构体
	std::deque<int> hoffmanCode;	//编码后序列
	int ch;					//字符
};

typedef struct Hnode
{ //霍夫曼树节点
	int weight;
	Hnode *lchild;
	Hnode *rchild;
	Hnode *parent;
	Hnode() :weight(0), lchild(NULL), rchild(NULL), parent(NULL) {}
};

struct cmp
{ //节点比较规则
	bool operator()(const Hnode *a, const Hnode *b)
	{
		return a->weight > b->weight;
	}
};

void getHuffmanCode(std::priority_queue<Hnode *, std::vector<Hnode *>, cmp> &, ChCode *code, Hnode *node, int n);
void getHuffmanTable(ChCode *code, int n);

#endif
