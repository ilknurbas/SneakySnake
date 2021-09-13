/* Include Files */
#include "SneakySnake.h"
#include "stdio.h"

/* Function Definitions */

/*
 * Arguments    : int ReadLength
 *                const char RefSeq[ReadLength]
 *                const char ReadSeq[ReadLength]
 *                int EditThreshold
 *                int KmerSize (t)
 *                int DebugMode
 *				  int IterationNo (y)
 * Return Type  : int Accepted
 */

 struct rowCharacteristics
 {
     char rowName[20];
     int rowNumber;
 } row;

int SneakySnake(int ReadLength, char * RefSeq, char * ReadSeq, int EditThreshold, int KmerSize, int DebugMode, int IterationNo, int * Checkpoint, int *rowNumber, char **cigarSSWithoutMatch)
{
	int Accepted=1;
	int n;
	int e;
	int K;
	int index=0;
	int Edits=0;
	int count = 0;
	int GlobalCount=0;
	
	int KmerStart=0;
	int KmerEnd=0;
	
	int roundsNo=0;
    row.rowNumber=0;


	if (DebugMode){
		printf("%2d - ", 0);
		for (n = 0; n < (ReadLength) ; n++) {
			if (n%10==0)
				printf("_%2d       ",n);
		}
		printf("\n");
		// Main Diagonal
		////////////////////////////////////////////
		printf("%2d - ", 0);
		for (n = 0; n < (ReadLength) ; n++) {
			if (ReadSeq[n]!= RefSeq[n]) {
				printf("1");
			}
			else if (ReadSeq[n]== RefSeq[n]) {
				printf("0");
			}
		}
		printf("\n");
		////////////////////////////////////////////
		//  Upper & Lower Diagonals
		////////////////////////////////////////////
		for (e = 1; e <= EditThreshold; e++) {
			count=0;
			///////////////////////////////////////////////////
			//  Shift Read to Right-hand side (Upper Diagonals: Deletion)
			printf("%2d - ",e);
			for (n = 0; n < (ReadLength) ; n++) {
				if (n<e) 
					printf("1");
				else if (ReadSeq[n-e]!= RefSeq[n])
					printf("1");
				else if (ReadSeq[n-e]== RefSeq[n]) {
					printf("0");
				}
			}
			printf("\n");
			///////////////////////////////////////////////////
			//  Shift Read to Left-hand side (Lower Diagonals: Insertion)
			printf("%2d - ",e+EditThreshold);
			for (n =0; n < (ReadLength) ; n++) {
				if (n>ReadLength-e-1) 
					printf("1");
				else if (ReadSeq[n+e]!= RefSeq[n])
					printf("1");
				else if (ReadSeq[n+e]== RefSeq[n]) {
					printf("0");
				}					
			}
			printf("\n");
			
		}
	}
	else{

	int Runtime = 1;

	    if (!Runtime) {
            // Go through each Kmer
            for(K=0; K < (ReadLength/KmerSize); K++){
                KmerStart= K*KmerSize;
                if (K < (ReadLength/KmerSize)-1)
                    KmerEnd = (K+1)*KmerSize;
                else
                    KmerEnd = ReadLength;

                index=(KmerStart);
                roundsNo=1;

                while (index<KmerEnd) {
                    fprintf(stderr, "[[SS::%s]] roundsNo: %d\n", __func__,roundsNo);
                    GlobalCount=0;

                    // Main Diagonal
                    ////////////////////////////////////////////
                    for (n = (index); n < (KmerEnd) ; n++) {
                        if (ReadSeq[n]!= RefSeq[n]) {
                            goto EXIT1;
                        }
                        else if (ReadSeq[n]== RefSeq[n]) {
                            GlobalCount=GlobalCount+1;
                        }
                    }
                    EXIT1:
                    if (GlobalCount == (KmerEnd-KmerStart)) {

                        fprintf(stderr, "[SS::%s] KmerEnd: %d\n", __func__,KmerEnd);
                        fprintf(stderr, "[SS::%s] KmerStart: %d\n", __func__,KmerStart);
                        fprintf(stderr, "[SS::%s] GlobalCount: %d\n", __func__,GlobalCount);
                        fprintf(stderr, "[SS::%s] Inside loop \n", __func__);

                        memcpy(row.rowName, "Main", 5);
                        row.rowName [4] = '\0';
                        row.rowNumber= 0;
                        *Checkpoint = GlobalCount;
                        goto LOOP_;
                    }

                    ////////////////////////////////////////////
                    //  Upper & Lower Diagonals
                    ////////////////////////////////////////////
                    for (e = 1; e <= EditThreshold; e++) {
                        count=0;
                        ///////////////////////////////////////////////////
                        //  Shift Read to Right-hand side (Upper Diagonals: Deletion)
                        for (n = (index); n < (KmerEnd) ; n++) {
                            if (n<e)
                                goto EXIT2; // fill the shifted chars with Ones
                            else if (ReadSeq[n-e]!= RefSeq[n])
                                goto EXIT2;
                            else if (ReadSeq[n-e]== RefSeq[n]) {
                                count=count+1;
                            }
                        }
                        EXIT2:
                        if (count>GlobalCount) {
                            GlobalCount = count;
                            row.rowNumber= e;
                            memcpy(row.rowName, "Deletion/Upper", 15);
                            row.rowName [14] = '\0';
                        }
                        if (count == (KmerEnd-KmerStart)){
                            goto LOOP;
                        }

                        count=0;
                        ///////////////////////////////////////////////////
                        //  Shift Read to Left-hand side (Lower Diagonals: Insertion)
                        for (n = (index); n < (KmerEnd) ; n++) {
                            if (n>ReadLength-e-1)
                                goto EXIT3;
                            else if (ReadSeq[n+e]!= RefSeq[n])
                                goto EXIT3;
                            else if (ReadSeq[n+e]== RefSeq[n]) {
                                count=count+1;
                            }
                        }
                        EXIT3:
                        if (count>GlobalCount) {
                            GlobalCount = count;
                            row.rowNumber= e + EditThreshold;
                            memcpy(row.rowName, "Insertion/Lower", 16);
                            row.rowName [15] = '\0';
                        }
                        if (count == (KmerEnd-KmerStart) ){
                            goto LOOP;
                        }

                    }

                    fprintf(stderr, "[SS::%s] Index before: %d\n", __func__,index);
                    index = index+GlobalCount; // we add one here to skip the error that causes the segmentation
                    if (index<(KmerEnd)) {
                        Edits=Edits+1;
                        index=index+1;
                    }

                    *Checkpoint = index-1;

                    fprintf(stderr, "[SS::%s] GlobalCount: %d\n", __func__,GlobalCount);
                    fprintf(stderr, "[SS::%s] Index after: %d\n", __func__,index);
                    fprintf(stderr, "[SS::%s] Checkpoint: %d\n", __func__,*Checkpoint);
                    fprintf(stderr, "[SS::%s] Edits: %d\n", __func__,Edits);
                    fprintf(stderr, "[SS::%s] Row Number: %d\n", __func__,row.rowNumber);

                    if (row.rowNumber == 0){
                          memcpy(row.rowName, "Main", 5);
                          row.rowName [4] = '\0';
                    }

                    if (roundsNo>IterationNo){
                      goto LOOP;
                    }
                    roundsNo=roundsNo+1;
                    if (Edits > EditThreshold)
                        return 0;


                    ////////////////////////////////////////////
                    // END of Building the Hamming masks
                    /////////////////////////////////////////////////
                    /////////////////////////////////////////////////
                    //// if not sure about the number of matches, try finding the best PATH within a Kmer

                }
                LOOP_:
                if (Edits > EditThreshold){
                    fprintf(stderr, "[SS::%s] Exceed!\n", __func__);
                    return 0;
                }

                *cigarSSWithoutMatch = NULL;
                fprintf(stderr, "[SS::%s] cigarSSWithoutMatch: %s\n", __func__,*cigarSSWithoutMatch);

                return 1;

                LOOP:
                /*
                GlobalCount=GlobalCount + Edits;
                printf("Global Count: %d, Edits: %d\n",GlobalCount,Edits);
                */

                if (Edits > EditThreshold){
                    fprintf(stderr, "[SS::%s] Exceed!\n", __func__);
                    return 0;
                }

                *rowNumber = row.rowNumber;
                int counter = row.rowNumber;
                int noOfDigits = 0;
                if(counter == 0){
                    noOfDigits =1;
                }else{
                    while (counter!=0) {
                        counter = counter/10;
                        noOfDigits++;
                    }
                }

                fprintf(stderr, "[SS::%s] Row Name: %s\n", __func__,row.rowName);
                char * returnValue = (char *)malloc( sizeof(char) * (noOfDigits));

                if (strcmp(row.rowName,"Main") == 0 ) {
                    sprintf(returnValue, "%d", 1);
                    returnValue = (char *) realloc(returnValue, sizeof(char) * (strlen(returnValue) + 2) );
                    memcpy (returnValue + strlen (returnValue), "X", 2 );
                    returnValue [ strlen(returnValue)] = '\0';

                }else if (strcmp(row.rowName,"Insertion/Lower") == 0 ) {
                    sprintf(returnValue, "%d", row.rowNumber);
                    returnValue = (char *) realloc(returnValue, sizeof(char) * (strlen(returnValue) + 2));
                    memcpy (returnValue + strlen (returnValue), "I", 2 );
                    returnValue [ strlen(returnValue)] = '\0';

                }else if (strcmp(row.rowName,"Deletion/Upper") == 0  ){
                    sprintf(returnValue, "%d", row.rowNumber);
                    returnValue = (char *) realloc(returnValue, sizeof(char) * (strlen(returnValue) + 2));
                    memcpy (returnValue + strlen (returnValue), "D", 2 );
                    returnValue [ strlen(returnValue)] = '\0';
                }

                // comes from parameter
                *cigarSSWithoutMatch = (char *)malloc( (strlen(returnValue) + 1) * sizeof(char));
                memcpy (*cigarSSWithoutMatch , returnValue, strlen(returnValue)  );
                (*cigarSSWithoutMatch)[strlen(returnValue)] = '\0';
                free(returnValue);

                fprintf(stderr, "[SS::%s] *cigarSSWithoutMatch: %s\n", __func__, *cigarSSWithoutMatch);

            }

		} else { // no output
            // Go through each Kmer
            for(K=0; K < (ReadLength/KmerSize); K++){
                KmerStart= K*KmerSize;
                if (K < (ReadLength/KmerSize)-1)
                    KmerEnd = (K+1)*KmerSize;
                else
                    KmerEnd = ReadLength;

                index=(KmerStart);
                roundsNo=1;

                while (index<KmerEnd) {
                    GlobalCount=0;

                    // Main Diagonal
                    ////////////////////////////////////////////
                    for (n = (index); n < (KmerEnd) ; n++) {
                        if (ReadSeq[n]!= RefSeq[n]) {
                            goto ExIT1;
                        }
                        else if (ReadSeq[n]== RefSeq[n]) {
                            GlobalCount=GlobalCount+1;
                        }
                    }
                    ExIT1:
                    if (GlobalCount == (KmerEnd-KmerStart)) {

                        memcpy(row.rowName, "Main", 5);
                        row.rowName [4] = '\0';
                        row.rowNumber= 0;
                        *Checkpoint = GlobalCount;
                        goto LoOP_;
                    }

                    ////////////////////////////////////////////
                    //  Upper & Lower Diagonals
                    ////////////////////////////////////////////
                    for (e = 1; e <= EditThreshold; e++) {
                        count=0;
                        ///////////////////////////////////////////////////
                        //  Shift Read to Right-hand side (Upper Diagonals: Deletion)
                        for (n = (index); n < (KmerEnd) ; n++) {
                            if (n<e)
                                goto ExIT2; // fill the shifted chars with Ones
                            else if (ReadSeq[n-e]!= RefSeq[n])
                                goto ExIT2;
                            else if (ReadSeq[n-e]== RefSeq[n]) {
                                count=count+1;
                            }
                        }
                        ExIT2:
                        if (count>GlobalCount) {
                            GlobalCount = count;
                            row.rowNumber= e;
                            memcpy(row.rowName, "Deletion/Upper", 15);
                            row.rowName [14] = '\0';
                        }
                        if (count == (KmerEnd-KmerStart)){
                            goto LoOP;
                        }

                        count=0;
                        ///////////////////////////////////////////////////
                        //  Shift Read to Left-hand side (Lower Diagonals: Insertion)
                        for (n = (index); n < (KmerEnd) ; n++) {
                            if (n>ReadLength-e-1)
                                goto ExIT3;
                            else if (ReadSeq[n+e]!= RefSeq[n])
                                goto ExIT3;
                            else if (ReadSeq[n+e]== RefSeq[n]) {
                                count=count+1;
                            }
                        }
                        ExIT3:
                        if (count>GlobalCount) {
                            GlobalCount = count;
                            row.rowNumber= e + EditThreshold;
                            memcpy(row.rowName, "Insertion/Lower", 16);
                            row.rowName [15] = '\0';
                        }
                        if (count == (KmerEnd-KmerStart) ){
                            goto LoOP;
                        }

                    }

                    index = index+GlobalCount; // we add one here to skip the error that causes the segmentation
                    if (index<(KmerEnd)) {
                        Edits=Edits+1;
                        index=index+1;
                    }

                    *Checkpoint = index-1;

                    if (row.rowNumber == 0){
                          memcpy(row.rowName, "Main", 5);
                          row.rowName [4] = '\0';
                    }

                    if (roundsNo>IterationNo){
                      goto LoOP;
                    }
                    roundsNo=roundsNo+1;
                    if (Edits > EditThreshold)
                        return 0;

                }
                LoOP_:
                if (Edits > EditThreshold){
                    return 0;
                }

                *cigarSSWithoutMatch = NULL;
                return 1;

                LoOP:
                /*
                GlobalCount=GlobalCount + Edits;
                printf("Global Count: %d, Edits: %d\n",GlobalCount,Edits);
                */

                if (Edits > EditThreshold){
                    return 0;
                }

                *rowNumber = row.rowNumber;
                int counter = row.rowNumber;
                int noOfDigits = 0;
                if(counter == 0){
                    noOfDigits =1;
                }else{
                    while (counter!=0) {
                        counter = counter/10;
                        noOfDigits++;
                    }
                }
                char * returnValue = (char *)malloc( sizeof(char) * (noOfDigits));

                if (strcmp(row.rowName,"Main") == 0 ) {
                    sprintf(returnValue, "%d", 1);
                    returnValue = (char *) realloc(returnValue, sizeof(char) * (strlen(returnValue) + 2) );
                    memcpy (returnValue + strlen (returnValue), "X", 2 );
                    returnValue [ strlen(returnValue)] = '\0';

                }else if (strcmp(row.rowName,"Insertion/Lower") == 0 ) {
                    sprintf(returnValue, "%d", row.rowNumber);
                    returnValue = (char *) realloc(returnValue, sizeof(char) * (strlen(returnValue) + 2));
                    memcpy (returnValue + strlen (returnValue), "I", 2 );
                    returnValue [ strlen(returnValue)] = '\0';

                }else if (strcmp(row.rowName,"Deletion/Upper") == 0  ){
                    sprintf(returnValue, "%d", row.rowNumber);
                    returnValue = (char *) realloc(returnValue, sizeof(char) * (strlen(returnValue) + 2));
                    memcpy (returnValue + strlen (returnValue), "D", 2 );
                    returnValue [strlen(returnValue)] = '\0';
                }

                // comes from parameter
                *cigarSSWithoutMatch = (char *)malloc( (strlen(returnValue) + 1) * sizeof(char));
                memcpy (*cigarSSWithoutMatch , returnValue, strlen(returnValue)  );
                (*cigarSSWithoutMatch)[strlen(returnValue)] = '\0';
                free(returnValue);

            }

		}
	}

	return Accepted;
}
