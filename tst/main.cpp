#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Database.h"
#include "DTW.h"
#include "DCRC.h"
#include "KEOGH.h"
#include "DCRC_Tree.h"
#include "RTree.h"
#include "time.h"
#include <malloc.h>
#include <crtdbg.h>
#include "debugleak.h"
#include "math.h"

#define LAMBDA 1
#define EPSILON 0.1
#define TOTAL_DATA 5504
#define SAMPLE_COUNT 100000
#define COPY_NUM (SAMPLE_COUNT/TOTAL_DATA)

#define TABLE2_FLAG 0
#define TABLE3_FLAG 0

#define CYLND_CURVE_TREE 2
#define TIGHTNESS_PRUNNING 0
#define LONG_SHORT 0
#define LAMBDA_TOLERANCE 1
#define DIM_FLAG 0
#define TEST_FLAG 0

#define TEX_DATA_BEGIN_1 "\\addplot[green, mark=*] plot coordinates{"
#define TEX_DATA_BEGIN_2 "\\addplot[blue, mark=square*] plot coordinates{"
#define TEX_DATA_BEGIN_3 "\\addplot[orange, mark=o] plot coordinates{"
#define TEX_DATA_BEGIN_4 "\\addplot[purple, mark=+] plot coordinates{"
#define TEX_DATE_END "};"
#define TEST_COUNT 10

using namespace std;

int DIM_LONG    = 20;
int DIM_SHORT   = 17;

#define CLASS_TEST_NUM 1

Database db;

void addTexBegin(int i) {
    char beginStr[100];
    char command[200];
    switch (i) {
    case 1:
        system("del p*.txt");
        sprintf_s(beginStr, TEX_DATA_BEGIN_1);
        break;
    case 2:
        sprintf_s(beginStr, TEX_DATA_BEGIN_2);
        break;
    case 3:
        sprintf_s(beginStr, TEX_DATA_BEGIN_3);
        break;
    case 4:
        sprintf_s(beginStr, TEX_DATA_BEGIN_4);
        break;
    }
    for (int j = 1; j <= 4; j ++) {
        sprintf_s(command, "echo %s >> p%d.txt", beginStr, j);
        system(command);
    }
}
void addTexEnd() {
    char command[200];
    for (int j = 1; j <= 4; j ++) {
        sprintf_s(command, "echo %s >> p%d.txt", TEX_DATE_END, j);
        system(command);
    }
}

void addTexBegin(int i, int p) {
    char beginStr[100];
    char command[200];
    switch (i) {
    case 1:
        sprintf_s(command, "del p%d.txt", p);
        system(command);
        sprintf_s(beginStr, TEX_DATA_BEGIN_1);
        break;
    case 2:
        sprintf_s(beginStr, TEX_DATA_BEGIN_2);
        break;
    case 3:
        sprintf_s(beginStr, TEX_DATA_BEGIN_3);
        break;
    case 4:
        sprintf_s(beginStr, TEX_DATA_BEGIN_4);
        break;
    }
    sprintf_s(command, "echo %s >> p%d.txt", beginStr, p);
    system(command);
}

void addTexEnd(int p) {
    char command[200];
    sprintf_s(command, "echo %s >> p%d.txt", TEX_DATE_END, p);
    system(command);
}



float dtw_one2Set(TimeSeries *one, TimeSeries** tmsSet, int nSet, float lambda) {
    float ret = INFINITY;
    for (int i = 0; i < nSet; i ++) {
        TimeSeries* two = tmsSet[i];
        DTW dtw;
        float tmp = dtw.calculate(one, two, RTree::mode, lambda);
        if (tmp < ret) ret = tmp;


    }
    return ret;
}


void calcDCRC(TimeSeries **testData, int index, float lambda, float &error) {
    Dataset dSet = db.datasets[index];
    int nClass = dSet.classCount;
    float errSum = 0.0;
    int errCnt = 0;
   for (int j = -1; j <= nClass; j ++) {
        int nSet;
        TimeSeries** tSet = db.getInnerClass(index,j,DIM_LONG,&nSet);
        if (nSet == 0) continue;
        DCRC dcrc(DIM_LONG, DIM_SHORT);
        dcrc.calculate(tSet, nSet);
        for (int k = 0; k < TEST_COUNT; k ++) {
            float tmpError= dcrc.query(testData[k], lambda);
            errSum += tmpError;
            errCnt ++;
        }
        delete tSet;
    }
    error = errSum / errCnt;
 }

void calcKEOGH(TimeSeries **testData, int index, float lambda, float &error) {
    Dataset dSet = db.datasets[index];
    int nClass = dSet.classCount;
    float errSum = 0.0;
    int errCnt = 0;
    for (int j = -1; j <= nClass; j ++) {
        int nSet;
        TimeSeries** tSet = db.getInnerClass(index,j,DIM_LONG,&nSet);
        if (nSet == 0) continue;
        KEOGH keogh;
        keogh.calculate(tSet, nSet, lambda);
        for (int k = 0; k < TEST_COUNT; k ++) {
            float tmpError= keogh.query(testData[k], lambda);
            errSum += tmpError;
            errCnt ++;
        }
        delete tSet;
    }
    error = errSum / errCnt;
}

void calcDCRCTightness(TimeSeries **testData, int index, float lambda, float &error) {
    Dataset dSet = db.datasets[index];
    int nClass = dSet.classCount;
    float errSum = 0.0;
    int errCnt = 0;
   for (int j = -1; j <= nClass; j ++) {
        int nSet;
        TimeSeries** tSet = db.getInnerClass(index,j,DIM_LONG,&nSet);
        if (nSet == 0) continue;
        DCRC dcrc(DIM_LONG, DIM_SHORT);
        dcrc.calculate(tSet, nSet);
        for (int k = 0; k < TEST_COUNT; k ++) {
            float tmpError = dcrc.query(testData[k], lambda);
            float dtwError = dtw_one2Set(testData[k], tSet, nSet, lambda);
            errSum += tmpError/dtwError;
            errCnt ++;
        }
        delete tSet;
    }
    error = errSum / errCnt;
 }

void calcKEOGHTightness(TimeSeries **testData, int index, float lambda, float &error) {
    Dataset dSet = db.datasets[index];
    int nClass = dSet.classCount;
    float errSum = 0.0;
    int errCnt = 0;
    for (int j = -1; j <= nClass; j ++) {
        int nSet;
        TimeSeries** tSet = db.getInnerClass(index,j,DIM_LONG,&nSet);
        if (nSet == 0) continue;
        KEOGH keogh;
        keogh.calculate(tSet, nSet, lambda);
        for (int k = 0; k < TEST_COUNT; k ++) {
            float tmpError= keogh.query(testData[k], lambda);
            float dtwError = dtw_one2Set(testData[k], tSet, nSet, lambda);
            errSum += tmpError/dtwError;
            errCnt ++;
        }
        delete tSet;
    }
    error = errSum / errCnt;
}

void calcDCRCPrune(TimeSeries **testData, int index, float lambda, float tolerance, float &error) {
    Dataset dSet = db.datasets[index];
    int nClass = dSet.classCount;
    int errCnt = 0;
    for (int j = -1; j <= nClass; j ++) {
        int nSet;
        TimeSeries** tSet = db.getInnerClass(index,j,DIM_LONG,&nSet);
        if (nSet == 0) continue;
        DCRC dcrc(DIM_LONG, DIM_SHORT);
        dcrc.calculate(tSet, nSet);
        for (int k = 0; k < TEST_COUNT; k ++) {
            float tmpError = dcrc.query(testData[k], lambda);
            if (tmpError > tolerance) errCnt++;
        }
        delete tSet;
    }
    error = 1.0f * errCnt / (TEST_COUNT * nClass);
 }

void calcKEOGHPrune(TimeSeries **testData, int index, float lambda, float tolerance, float &error) {
    Dataset dSet = db.datasets[index];
    int nClass = dSet.classCount;
    int errCnt = 0;
    for (int j = -1; j <= nClass; j ++) {
        int nSet;
        TimeSeries** tSet = db.getInnerClass(index,j,DIM_LONG,&nSet);
        if (nSet == 0) continue;
        KEOGH keogh;
        keogh.calculate(tSet, nSet, lambda);
        for (int k = 0; k < TEST_COUNT; k ++) {
            float tmpError= keogh.query(testData[k], lambda);
            if (tmpError > tolerance) errCnt++;
        }
        delete tSet;
    }
    error = 1.0f * errCnt / (TEST_COUNT * nClass);
}



void calcDTW(TimeSeries **testData, int index, float lambda, float &error) {
    Dataset dSet = db.datasets[index];
    int nClass = dSet.classCount;
    float errSum = 0.0;
    int errCnt = 0;
    for (int j = -1; j <= nClass; j ++) {
        int nSet;
        TimeSeries** tSet = db.getInnerClass(index,j,DIM_LONG,&nSet);
        if (nSet == 0) continue;
        for (int k = 0; k < TEST_COUNT; k ++) {
            float tmpError= dtw_one2Set(testData[k], tSet, nSet, lambda);
            errSum += tmpError;
            errCnt ++;
        }
        delete tSet;
    }
    error = errSum / errCnt;
}
void calcCurve() {
    TimeSeries ** testData = db.getTestData(DIM_LONG);
    TimeSeries ** testDataList[7];
    char command[200];
    if (DIM_FLAG) {
        for (int i = 0; i <= 6; i ++) {
            testDataList[i] = db.getTestData(20 + 5 * i);
        }
        if (LONG_SHORT == 0) {
            for (int index = 0; index < 4; index++) {
                printf("data index = %d: \n", index);
                addTexBegin(1, index);
                for (int i = 0; i <= 6; i ++) {
                    testData = testDataList[i];
                    RTree::mode = 0;
                    float keoghError = 0;
                    DIM_LONG = 20 + 5 * i;
                    calcKEOGHTightness(testData, index, (float)LAMBDA, keoghError);
                    sprintf_s(command, "echo (%d,\t%0.3f) >> p%d.txt", DIM_LONG, keoghError, index);
                    system(command);
                    printf("(%d,\t%0.3f)\n", DIM_LONG, keoghError);
                }
                addTexEnd(index);

                printf("\n");
                addTexBegin(2, index);
                for (int i = 0; i <= 6; i ++) {
                    testData = testDataList[i];
                    RTree::mode = 1;
                    float cornerError = 0;
                    DIM_LONG = 20 + 5 * i;
                    calcKEOGHTightness(testData, index, (float)LAMBDA, cornerError);
                    sprintf_s(command, "echo (%d,\t%0.3f) >> p%d.txt", DIM_LONG, cornerError, index);
                    system(command);
                    printf("(%d,\t%0.3f)\n", DIM_LONG, cornerError);
                }
                addTexEnd(index);


                addTexBegin(3, index);
                printf("\n");
                for (int i = 0; i <= 6; i ++) {
                    testData = testDataList[i];
                    RTree::mode = 0;
                    DCRC_Tree::mode = 0;
                    float dcrcError = 0;
                    DIM_LONG = 20 + 5 * i;
                    calcDCRCTightness(testData, index, (float)LAMBDA, dcrcError);
                    sprintf_s(command, "echo (%d,\t%0.3f) >> p%d.txt", DIM_LONG, dcrcError, index);
                    system(command);
                    printf("(%d,\t%0.3f)\n", DIM_LONG, dcrcError);
                }
                addTexEnd(index);
            }
        }
        else {
            for (int index = 0; index < 4; index++) {
                printf("data index = %d: \n", index);
                addTexBegin(index + 1,0);
                for (int i = 0; i <= 6; i ++) {
                    DCRC_Tree::mode = 0;
                    float dcrcError = 0;
                    DIM_SHORT = 16 - 2 * i;
                    calcDCRCTightness(testData, index, (float)LAMBDA, dcrcError);
                    sprintf_s(command, "echo (%d,\t%0.3f) >> p0.txt", DIM_SHORT, dcrcError);
                    system(command);
                    printf("(%d,\t%0.3f)\n", DIM_SHORT, dcrcError);
                }
                addTexEnd(0);
            }
            for (int index = 0; index < 4; index++) {
                printf("data index = %d: \n", index + 4);
                addTexBegin(index + 1,1);
                for (int i = 0; i <= 6; i ++) {
                    DCRC_Tree::mode = 0;
                    float dcrcError = 0;
                    DIM_SHORT = 16 - 2 * i;
                    calcDCRCTightness(testData, index + 4, (float)LAMBDA, dcrcError);
                    sprintf_s(command, "echo (%d,\t%0.3f) >> p1.txt", DIM_SHORT, dcrcError);
                    system(command);
                    printf("(%d,\t%0.3f)\n", DIM_SHORT, dcrcError);
                }
                addTexEnd(1);
            }
        }
    }
    else if (TIGHTNESS_PRUNNING == 0) {
        for (int index = 0; index < 4; index++) {
            addTexBegin(1, index);
            printf("\nindex=%d\n", index);
            for (int i = 1; i <= 10; i ++) {
                float lambda = 0.1f * i;
                RTree::mode = 0;
                float keoghError = 0;
                calcKEOGHTightness(testData, index, lambda, keoghError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", lambda, keoghError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", lambda, keoghError);
            }
            addTexEnd(index);
            addTexBegin(2, index);
            printf("\n");
            for (int i = 1; i <= 10; i ++) {
                float lambda = 0.1f * i;
                RTree::mode = 1;
                float cornerError = 0;
                calcKEOGHTightness(testData, index, lambda, cornerError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", lambda, cornerError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", lambda, cornerError);
            }
            addTexEnd(index);

            addTexBegin(3, index);
            printf("\n");
            for (int i = 1; i <= 10; i ++) {
                float lambda = 0.1f * i;
                DCRC_Tree::mode = 0;
                float dcrcError = 0;
                calcDCRCTightness(testData, index, lambda, dcrcError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", lambda, dcrcError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", lambda, dcrcError);
            }
            addTexEnd(index);
        }
    }
    else if (LAMBDA_TOLERANCE == 0) {
        for (int index = 0; index < 4; index++) {
            printf("\nindex=%d\n", index);

            addTexBegin(1, index);
            for (int i = 1; i <= 10; i ++) {
                float lambda = 0.1f * i;
                RTree::mode = 0;
                float keoghError = 0;
                calcKEOGHPrune(testData, index, lambda, 0.1f, keoghError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", lambda, keoghError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", lambda, keoghError);
            }
            addTexEnd(index);

            printf("\n");
            addTexBegin(2, index);
            for (int i = 1; i <= 10; i ++) {
                float lambda = 0.1f * i;
                RTree::mode = 1;
                float cornerError = 0;
                calcKEOGHPrune(testData, index, lambda, 0.1f, cornerError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", lambda, cornerError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", lambda, cornerError);
            }
            addTexEnd(index);

            printf("\n");
            addTexBegin(3, index);
            for (int i = 1; i <= 10; i ++) {
                float lambda = 0.1f * i;
                DCRC_Tree::mode = 0;
                float dcrcError = 0;
                calcDCRCPrune(testData, index, lambda, 0.1f, dcrcError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", lambda, dcrcError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", lambda, dcrcError);
            }
            addTexEnd(index);
        }
    }
    else {
        for (int index = 0; index < 4; index++) {
            printf("\nindex=%d:\n", index);
            addTexBegin(1, index);
            for (int i = 1; i <= 10; i ++) {
                float tolerance = 0.1f * i;
                RTree::mode = 0;
                float keoghError = 0;
                calcKEOGHPrune(testData, index, (float)LAMBDA, tolerance, keoghError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", tolerance, keoghError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", tolerance, keoghError);
            }
            addTexEnd(index);

            printf("\n");
            addTexBegin(2, index);
            for (int i = 1; i <= 10; i ++) {
                float tolerance = 0.1f * i;
                RTree::mode = 1;
                float cornerError = 0;
                calcKEOGHPrune(testData, index, (float)LAMBDA, tolerance, cornerError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", tolerance, cornerError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", tolerance, cornerError);
            }
            addTexEnd(index);

            addTexBegin(3, index);
            printf("\n");
            for (int i = 1; i <= 10; i ++) {
                float tolerance = 0.1f * i;
                DCRC_Tree::mode = 0;
                float dcrcError = 0;
                calcDCRCPrune(testData, index, (float)LAMBDA, tolerance, dcrcError);
                sprintf_s(command, "echo (%0.1f,\t%0.3f) >> p%d.txt", tolerance, dcrcError, index);
                system(command);
                printf("(%0.1f,\t%0.3f)\n", tolerance, dcrcError);
            }
            addTexEnd(index);
        }
    }

}
void calcCylindar(float error) {
    TimeSeries ** testData = db.getTestData(DIM_LONG);
    for (int i = 0; i < 20; i ++) {
        if (TABLE2_FLAG) {
            //table2
            RTree::mode = 0;
            DCRC_Tree::mode = 0;
            float dcrcError1 = 0;
            calcDCRC(testData, i, 0.2f, dcrcError1);

            float dtwError1 = 0;
            calcDTW(testData, i, 0.2f, dtwError1);

            RTree::mode = 0;
            DCRC_Tree::mode = 0;
            float dcrcError2 = 0;
            calcDCRC(testData, i, 0.6f, dcrcError2);

            float dtwError2 = 0;
            calcDTW(testData, i, 0.6f, dtwError2);

            RTree::mode = 0;
            DCRC_Tree::mode = 0;
            float dcrcError3 = 0;
            calcDCRC(testData, i, 1.0f, dcrcError3);

            float dtwError3 = 0;
            calcDTW(testData, i, 1.0f, dtwError3);

            printf("%s & %0.3f/%0.3f & %0.3f/%0.3f & %0.3f/%0.3f \\\\",
               db.datasets[i].dataName,
               dcrcError1,
               dtwError1,
               dcrcError2,
               dtwError2,
               dcrcError3,
               dtwError3
            );
        }
        else if (TABLE3_FLAG) {
            //table3
            float dtwError = 0;
            calcDTW(testData, i, 0.2f, dtwError);

            RTree::mode = 0;
            DCRC_Tree::mode = 0;
            float dcrcError = 0;
            calcDCRC(testData, i, 0.2f, dcrcError);

            RTree::mode = 0;
            float keoghError = 0;
            calcKEOGH(testData, i, 0.2f, keoghError);

            RTree::mode = 1;
            float cornerError = 0;
            calcKEOGH(testData, i, 0.2f, cornerError);

            printf("%s & %0.3f & %0.3f & %0.3f & %0.3f\\\\",
                db.datasets[i].dataName,
                dtwError,
                dcrcError,
                keoghError,
                cornerError
            );
        }
        else if (TIGHTNESS_PRUNNING == 0) {
            //tightness
            RTree::mode = 0;
            DCRC_Tree::mode = 0;
            float dcrcError = 0;
            calcDCRCTightness(testData, i, (float)LAMBDA, dcrcError);

            RTree::mode = 0;
            float keoghError = 0;
            calcKEOGHTightness(testData, i, (float)LAMBDA, keoghError);

            RTree::mode = 1;
            float cornerError = 0;
            calcKEOGHTightness(testData, i, (float)LAMBDA, cornerError);

/*            printf("%d\t%0.3f\t%0.3f\t%0.3f",
                i + 1,
                keoghError,
                cornerError,
                dcrcError
            );
*/
            printf("%d,%0.3f,%0.3f,%0.3f",
                i + 1,
                keoghError,
                cornerError,
                dcrcError
            );
        }
        else if (TIGHTNESS_PRUNNING == 1) {
            // prunning power
            RTree::mode = 0;
            DCRC_Tree::mode = 0;
            float dcrcError = 0;
            calcDCRCPrune(testData, i, (float)LAMBDA, (float)EPSILON, dcrcError);

            RTree::mode = 0;
            float keoghError = 0;
            calcKEOGHPrune(testData, i, (float)LAMBDA, (float)EPSILON, keoghError);

            RTree::mode = 1;
            float cornerError = 0;
            calcKEOGHPrune(testData, i, (float)LAMBDA, (float)EPSILON, cornerError);


/*            printf("%d\t%0.3f\t%0.3f\t%0.3f",
                i + 1,
                keoghError,
                cornerError,
                dcrcError
            );
*/
            printf("%d,%0.3f,%0.3f,%0.3f",
                i + 1,
                keoghError,
                cornerError,
                dcrcError
            );
        }
        printf("\n");
    }
    delete testData;
}

void queryOnDTree(DCRC_Tree * dtree) {
    clock_t t1 = clock();
    int n1, n2;
    int n3 = 0;
    for (int i = 0; i < 20; i ++) {
        for (int j = 0; j < CLASS_TEST_NUM; j ++) {
            TimeSeries *tms = &db.datasets[i].testData[j];
            //printf("query_%d_%d:\n", i, j);
            dtree -> query(tms, (float)LAMBDA, (float)EPSILON, n1, n2);

            DTW dtw;
            for (int k = 0; k < n2; k ++) {
                float tmp = dtw.calculate(tms, dtree -> resultList[k], DCRC_Tree::mode, (float)LAMBDA);
                if (tmp <= EPSILON) n3 ++;
            }
        }
    }
    clock_t t2 = clock();
    if (n3 > 0) printf("n3 = %d\n", n3);
    printf("d = %d, s = %d\n", DIM_LONG, t2 - t1);
    char command[100];
    sprintf_s(command, "echo (%d,%8.3f) >> p0.txt", DIM_LONG, (t2-t1)/(20.0 * CLASS_TEST_NUM));
    system(command);
}

double getUCRsuite(TimeSeries* a, Range* r) {
	double res = 0;
	for (int i = 0; i < MAX_DIM; i++) {
		if (a->data[i] > r[i].rMax) {
			double diff = (a->data[i] - r[i].rMax);
			res += diff * diff;
		}
		else if (a->data[i] < r[i].rMin){
			double diff = (r[i].rMin - a->data[i]);
			res += diff * diff;
		}
	}
	return sqrt(res);
}

/*************************************************************/
// There is a serious mistake in UCR suite !!!
// The algorithm is based on the following theoreom:
// DTW(Q(1:k),C(1:k)) + LB_Keogh(Q(k+1:n),C(k+1:n)) is a lower bound of DTW(Q(1:n),C(1:n))
// is not satisfied !!!
// A simple example: Q=1,1,2  C=1,2,2  n=3  k=2
/*************************************************************/
void searchByUCRsuite(TimeSeries *qts) {
	Range ranges[MAX_DIM];
	for (int i = 0; i < MAX_DIM; i++) {
		ranges[i].rMin = INFINITY;
		ranges[i].rMax = -INFINITY;
	}
	for (int i = 0; i < MAX_DIM; i++) {
		for (int j = 0; j < MAX_DIM; j++) {
			if ((abs(i - j) < LAMBDA*MAX_DIM)) {
				if (qts->data[i] > ranges[j].rMax) {
					ranges[j].rMax = qts->data[i];
				}
				else if (qts->data[i] < ranges[j].rMin) {
					ranges[j].rMin = qts->data[i];
				}
			}
		}
	}

	printf("\n%d results:", DCRC_Tree::resNodeNum);
	int sumSize = 0;
	for (int i = 0; i < DCRC_Tree::resNodeNum; i++) {
		printf("%d ", DCRC_Tree::resIndices[i]);
		char fileName[100] = { 0 };
		sprintf_s(fileName, "%s\\cluster\\%d.dat", DB_PATH, DCRC_Tree::resIndices[i]);
		FILE* fp = fopen(fileName, "rb");
		if (fp) {
			fseek(fp, 0, SEEK_END);
			int size = ftell(fp) / ( sizeof(float) * MAX_DIM );
			sumSize += size;
			TimeSeries ts;
			ts.length = MAX_DIM;
			fseek(fp, 0, 0);
			for (int j = 0; j < size; j ++) {
				fread(ts.data, sizeof(float), MAX_DIM, fp);
				getUCRsuite(&ts, ranges);
			}
			fclose(fp);
		}
	}
	printf("\n");
	printf("---------- %d, %d, %.1f%%\n", 
		sumSize, SAMPLE_COUNT, 100.0 * sumSize / SAMPLE_COUNT);
}
void queryOnDTree(DCRC_Tree *dtree, float lambda, float tolerance) {
    int n1, n2, n3;
    TimeSeries *tms = &db.datasets[0].trainData[0];
    dtree -> query(tms, 0.1f, INFINITY, n1, n2);
    //printf("%d\n", n1);
    int nodeNumber = n1;

    clock_t t1 = clock();
    float totalNodes = 0;
    float totalSeries = 0;
    float rightness = 0;
    for (int i = 0; i < 20; i ++) {
        for (int j = 0; j < CLASS_TEST_NUM; j ++) {
            TimeSeries *tms = &db.datasets[i].testData[j];
            //printf("i=%d, j= %d ", i, j);
            dtree -> query(tms, lambda, tolerance, n1, n2);
            //printf("n2 = %d\n", n2);
            totalNodes += n1;
            //totalSeries += n2;
			searchByUCRsuite(tms);
			/*
				DTW dtw;
				n3 = 0;
				for (int k = 0; k < n2; k ++) {
				float tmp = dtw.calculate(tms, dtree -> resultList[k], DCRC_Tree::mode, lambda);
				if (tmp <= tolerance) n3 ++;
				}
				rightness += n3;
			*/

        }
    }
    clock_t t2 = clock();

    char command[100];
    float param = ((LAMBDA_TOLERANCE == 0)? lambda : tolerance);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p1.txt", param, 100.0 * totalNodes/(20.0 * CLASS_TEST_NUM * nodeNumber));
    system(command);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p2.txt", param, 100.0 * totalSeries/(20.0 * CLASS_TEST_NUM * SAMPLE_COUNT));
    system(command);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p3.txt", param, 100.0 * rightness/totalSeries);
    system(command);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p4.txt", param, (t2 - t1)/(20.0 * CLASS_TEST_NUM * 60.0));
    system(command);

    printf("nodeNumber = %d\n", nodeNumber);
    printf("(%0.1f, %8.3f, %8.3f, %8.3f, %8.3f)\n",
        param,
        100.0 * totalNodes/ (DCRC_Tree::totalLeafNumber * 20 * CLASS_TEST_NUM),
        100.0 * totalSeries/ (SAMPLE_COUNT * 20 * CLASS_TEST_NUM),
        0.0f, //100.0 * rightness/totalSeries,
        (t2 - t1)/(20.0 * CLASS_TEST_NUM)
    );
}

TimeSeries* tmsList[600];
DCRC_Tree * buildDTree() {
    DCRC_Tree * dtree = new DCRC_Tree(DIM_LONG, DIM_SHORT);
    int count = 0;
    printf("\nloading:\n");
	srand(time(0));
	int loop = 0;
	int tempGroup = 200;
	while (loop++ < SAMPLE_COUNT) {
		int i = 13; //rand() % 20;
		/* if (loop < 10 * tempGroup) {
			for (int k = 0; k < 10; k++) {
				for (int p = 0; p < tempGroup; p++) {
					int j = rand() % (db.datasets[i].trainCount);
					TimeSeries *tms = &db.datasets[i].trainData[j];
					tms->index = count++;
					tmsList[p] = tms;
				}
				dtree->insertSeriesList(tmsList, tempGroup);
			}
		}
		else {
		*/
			int j = rand() % (db.datasets[i].trainCount);
			TimeSeries *tms = &db.datasets[i].trainData[j];
			dtree->insertSeries(tms);
		//}
		if (loop % 1000 == 0) 
			printf("\n%f%%\n", loop * 100.0 / SAMPLE_COUNT);
    }
    printf("\nnumber = %d\n\n", count);
    return dtree;
}
void calcDTree() {
    DCRC_Tree * dtree = buildDTree();
    
	queryOnDTree(dtree, (float)(LAMBDA), (float)EPSILON);
	/*
	for (int j = 1; j <= 10; j ++) {
	if (LAMBDA_TOLERANCE == 0)
	queryOnDTree(dtree, (float)(LAMBDA * j), (float)EPSILON);
	else
	queryOnDTree(dtree, (float)LAMBDA, (float)(j * EPSILON));
	}
	*/
    printf("\n\n");
}
void queryOnRTree(RTree * rtree) {
    clock_t t1 = clock();
    int n1, n2;
    int n3 = 0;
    for (int i = 0; i < 20; i ++) {
        for (int j = 0; j < CLASS_TEST_NUM; j ++) {
            TimeSeries *tms = &db.datasets[i].testData[j];
            //printf("i = %d, j = %d ", i, j);
            rtree -> query(tms, (float)LAMBDA, (float)EPSILON, n1, n2);
            //printf("n2 = %d\n", n2);
            /*
			DTW dtw;
			for (int k = 0; k < n2; k ++) {
			float tmp = dtw.calculate(tms, rtree -> resultList[k], RTree::mode, (float)LAMBDA);
			if (tmp <= EPSILON) n3 ++;
			}
			*/
        }
    }
    clock_t t2 = clock();
    if (n3 > 0) printf("n3 = %d\n", n3);
    printf("d = %d, s = %d\n", DIM_LONG, t2 - t1);
    char command[100];
    sprintf_s(command, "echo (%d,%8.3f) >> p0.txt", DIM_LONG, (t2-t1)/(20.0 * CLASS_TEST_NUM));
    system(command);
}

/*
TimeSeries * samples[TOTAL_DATA * COPY_NUM];
void buildSamples() {
int p = 0;
for (int i = 0; i < 20; i++) {
for (int k = 0; k < COPY_NUM; k++) {
for (int j = 0; j < db.datasets[i].trainCount; j++) {
TimeSeries *tms = &db.datasets[i].trainData[j];
samples[p++] = tms;
}
}
}
}
void queryLinear(float lambda, float epsilon) {
Range rngs[MAX_DIM];
clock_t t1 = clock();
DTW dtw;
KEOGH keogh;
int totalSeries = 0;
int rightness = 0;
for (int i = 0; i < 20; i++) {
for (int j = 0; j < CLASS_TEST_NUM; j++) {
TimeSeries *tms = &db.datasets[i].testData[j];
keogh.getRangesFromSeries(tms, lambda, rngs);
//printf("i = %d, j = %d ", i, j);
for (int k = 0; k < TOTAL_DATA * COPY_NUM; k++) {
if (keogh.calculate(samples[k], rngs) > epsilon) continue;
totalSeries++;
float tmp = dtw.calculate(tms, samples[k], RTree::mode, lambda);
if (tmp <= epsilon) rightness++;
}
}
}
clock_t t2 = clock();

char command[100];
float param = ((LAMBDA_TOLERANCE == 0) ? lambda : epsilon);
sprintf_s(command, "echo (%4.1f,%8.3f) >> p2.txt", param, 100.0 * totalSeries / (20.0 * CLASS_TEST_NUM * TOTAL_DATA * COPY_NUM));
system(command);
sprintf_s(command, "echo (%4.1f,%8.3f) >> p3.txt", param, 100.0 * rightness / totalSeries);
system(command);
sprintf_s(command, "echo (%4.1f,%8.3f) >> p4.txt", param, (t2 - t1) / (20.0 * CLASS_TEST_NUM * 60.0));
system(command);

printf("linear results: (%d, %8.3f)\n", rightness, (t2 - t1) / (20.0 * CLASS_TEST_NUM));
}
*/


void queryOnRTree(RTree * rtree, float lambda, float tolerance ) {
    int n1, n2, n3;

    TimeSeries *tms = &db.datasets[0].trainData[0];
    rtree -> query(tms, 0.1f, INFINITY, n1, n2);
    //printf("%d\n", n1);
    int nodeNumber = n1;

    clock_t t1 = clock();
    float totalNodes = 0;
    float totalSeries = 0;
    float rightness = 0;
    for (int i = 0; i < 20; i ++) {
        for (int j = 0; j < CLASS_TEST_NUM; j ++) {
            TimeSeries *tms = &db.datasets[i].testData[j];
            //printf("i = %d, j = %d ", i, j);
            rtree -> query(tms, lambda, tolerance, n1, n2);
            //printf("n2 = %d\n", n2);
            totalNodes += n1;
            totalSeries += n2;

            DTW dtw;
            n3 = 0;
            for (int k = 0; k < n2; k ++) {
                float tmp = dtw.calculate(tms, rtree -> resultList[k], RTree::mode, lambda);
                if (tmp <= tolerance) n3 ++;
            }
            rightness += n3;
        }
    }
    clock_t t2 = clock();

    char command[100];
    float param = ((LAMBDA_TOLERANCE == 0)? lambda : tolerance);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p1.txt", param, 100.0 * totalNodes/(20.0 * CLASS_TEST_NUM * nodeNumber));
    system(command);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p2.txt", param, 100.0 * totalSeries/(20.0 * CLASS_TEST_NUM * TOTAL_DATA * COPY_NUM));
    system(command);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p3.txt", param, 100.0 * rightness/totalSeries);
    system(command);
    sprintf_s(command, "echo (%4.1f,%8.3f) >> p4.txt", param, (t2 - t1)/(20.0 * CLASS_TEST_NUM * 60.0));
    system(command);

    printf("(%0.1f, %8.3f, %8.3f, %8.3f, %8.3f)\n",
        param,
        100.0 * totalNodes/(20.0 * CLASS_TEST_NUM * nodeNumber),
        100.0 * totalSeries/(20.0 * CLASS_TEST_NUM * TOTAL_DATA * COPY_NUM),
        100.0 * rightness/totalSeries,
        (t2 - t1)/ (20.0 * CLASS_TEST_NUM)
    );
}

RTree * buildRTree() {
    RTree * rtree = new RTree(DIM_LONG);
    int count = 0;
    printf("\nloading: ");
    for (int i = 0; i < 20; i ++) {
        printf("%d ", i);
        for (int k = 0; k < COPY_NUM; k ++) {
            for (int j = 0; j < db.datasets[i].trainCount; j ++) {
                TimeSeries *tms = &db.datasets[i].trainData[j];
                rtree -> insertSeries(tms, (float)LAMBDA);
                count ++;
            }
        }
    }
    printf("\nnumber = %d\n\n", count);
    return rtree;
}
void calcRTree() {
    RTree * rtree = buildRTree();
    for (int j = 1; j <= 10; j ++) {
      if (LAMBDA_TOLERANCE == 0)
        queryOnRTree(rtree, (float)(LAMBDA * j), (float)EPSILON);
      else
        queryOnRTree(rtree, (float)(LAMBDA), (float)(j * EPSILON));
    }
    printf("\n\n");
    //rtree -> print();
	delete rtree;
}

/*
	void calcLinear() {
	buildSamples();
	for (int j = 1; j <= 10; j++) {
	if (LAMBDA_TOLERANCE == 0)
	queryLinear((float)(LAMBDA * j), (float)EPSILON);
	else
	queryLinear((float)(LAMBDA), (float)(j * EPSILON));
	}
	printf("\n\n");
	}
*/


void compareTree() {
//	calcLinear();
	/*
	RTree * rtree = buildRTree();
	RTree::mode = 0;
	queryOnRTree(rtree, (float)(LAMBDA), (float)EPSILON);
	//rtree -> print();
	delete rtree;
	rtree = buildRTree();
	RTree::mode = 1;
	queryOnRTree(rtree, (float)(LAMBDA), (float)EPSILON);
	//rtree -> print();
	delete rtree;
	DCRC_Tree * dtree = buildDTree();
	queryOnDTree(dtree, (float)(LAMBDA), (float)EPSILON);
	printf("\n\n");
	delete dtree;

	buildSamples();
	RTree::mode = 0;
	queryLinear((float)(LAMBDA), (float)EPSILON);

	return;


    addTexBegin(1);
    RTree::mode = 0;
    calcRTree();
    addTexEnd();

    addTexBegin(2);
    RTree::mode = 1;
    calcRTree();
    addTexEnd();
	*/

	addTexBegin(3);
	RTree::mode = 0;
	DCRC_Tree::mode = 0;
	calcDTree();
	addTexEnd();

/*
	addTexBegin(4);
	RTree::mode = 0;
	calcLinear();
	addTexEnd();
*/

}

void calcTreeDim() {

    if (0) {
        addTexBegin(1, 0);
        for (int i = 0; i <= 6; i ++) {
            DIM_LONG = 20 + i * 5;
            RTree::mode = 1;
            RTree *rtree = buildRTree();
            queryOnRTree(rtree);
            delete rtree;
        }
        addTexEnd(0);

        addTexBegin(2, 0);
        for (int i = 0; i <= 6; i ++) {
            DIM_LONG = 20 + i * 5;
            RTree::mode = 1;
            RTree *rtree = buildRTree();
            queryOnRTree(rtree);
            delete rtree;
        }
        addTexEnd(0);
    }
    else {
        addTexBegin(3, 0);
        for (int i = 0; i <= 6; i ++) {
            DIM_LONG = 20 + i * 5;
            RTree::mode = 0;
            DCRC_Tree::mode = 0;
            DCRC_Tree *dtree = buildDTree();
            queryOnDTree(dtree);
            delete dtree;
        }
        addTexEnd(0);

    }
}

void testTree() {
    /*
    int n1, n2;
    TimeSeries *tms = db.datasets[0].trainData[0].dimTrans(DIM_LONG);
    //printf("%d\n", n1);
    RTree::mode = 0;
    RTree *rtree = buildRTree();
    rtree -> query(tms, LAMBDA, INFINITY, n1, n2);
    printf("number of non-leaf nodes of R-tree for LB_Keogh is %d\n", n1);

    RTree::mode = 1;
    rtree = buildRTree();
    rtree -> query(tms, LAMBDA, INFINITY, n1, n2);
    printf("number of non-leaf nodes of R-tree for LB_Corner is %d\n", n1);

    RTree::mode = 0;
    DCRC_Tree::mode = 0;
    //calcDTree();
    DCRC_Tree * dtree = buildDTree();
    dtree -> query(tms, LAMBDA, INFINITY, n1, n2);
    printf("number of non-leaf nodes of DCRC-tree is %d\n", n1);
    queryOnDTree(dtree, LAMBDA, EPSILON);
    printf("\n\n");
    */


    RTree::mode = 0;
    RTree *rtree = buildRTree();
    queryOnRTree(rtree, (float)LAMBDA, (float)EPSILON);



    RTree::mode = 1;
    rtree = buildRTree();
    queryOnRTree(rtree, (float)LAMBDA, (float)EPSILON);


    RTree::mode = 0;
    DCRC_Tree::mode = 0;
    DCRC_Tree * dtree = buildDTree();
    queryOnDTree(dtree, (float)LAMBDA, (float)EPSILON);

}

int main()
{
    cout << "Hello world!" << endl;
	//_CrtSetBreakAlloc(150); 
	db.setDim(DIM_LONG);
    db.read();

    cout << "Data is loaded successfully!" << endl;

    time_t t1 = time(0);
    
	switch (CYLND_CURVE_TREE) {
    case 0:
        calcCylindar(0.0);
        break;
    case 1:
        calcCurve();
        break;
    case 2:
        if (TEST_FLAG) {
            testTree();
        }
        else if (DIM_FLAG) {
            calcTreeDim();
       }
        else {
			printf("Compare trees: \n");
            compareTree();
        }
        break;
    }
    time_t t2 = time(0);
    printf("\n\nseconds = %d\n\n", (int)(t2 - t1));
	//   _CrtDumpMemoryLeaks();

   /*
    int nSet;
    TimeSeries** tSet1 = db.getInnerClass(5,1,DIM_LONG,&nSet);
    TimeSeries * ts1 = tSet1[0];
    TimeSeries * ts2 = tSet1[1];
    TimeSeries** tSet2 = db.getInnerClass(3,1,DIM_LONG,&nSet);
    TimeSeries * ts3 = tSet2[2];
    TimeSeries * ts4 = tSet2[3];
    DCRC dcrc1(DIM_LONG, DIM_SHORT);
    dcrc1.reset();
    dcrc1.addSeries(ts1);
    dcrc1.addSeries(ts2);

    DCRC dcrc2(DIM_LONG, DIM_SHORT);
    dcrc2.reset();
    dcrc2.addSeries(ts3);
    dcrc2.addSeries(ts4);
    //dcrc2.print();

    DCRC dcrc3(DIM_LONG, DIM_SHORT - 5);
    dcrc3.reset();
    dcrc3.addDCRC(&dcrc2);
    float t1 = dcrc3.testVolume(&dcrc1);
    printf("t1 = %f\n", t1);
    float t2 = dcrc3.testVolume(&dcrc2);
    printf("t2 = %f\n", t2);

    //dcrc3.print();
    */

/*
    //printf("------------------------------------------------------------\n");
    //calcDTW();
    printf("------------------------------------------------------------\n");
    calcDCRC(LAMBDA);
    printf("------------------------------------------------------------\n");
    calcKEOGH(0,LAMBDA);
    printf("------------------------------------------------------------\n");
    calcKEOGH(1,LAMBDA);
    printf("------------------------------------------------------------\n");
*/
    getchar();


    return 0;
}
