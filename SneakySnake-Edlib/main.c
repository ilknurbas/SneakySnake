/* Include Files */
#include "stdio.h" 
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <immintrin.h>
#include "SneakySnake.h"
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include "edlib.h"
#include "stdint.h"

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

int * Edlib_Edit ;
int * SS_Edlib_Edit ;

typedef struct {
	FILE *fp;
	int max_lines, buf_size, n_threads;
	int DebugMode,KmerSize, ReadLength,IterationNo, TotalAccepted, EditThreshold;
	char *buf;
} pipeline_t;

typedef struct {
	int n_lines;
	int DebugMode,KmerSize, ReadLength, IterationNo, EditThreshold;
	char **lines;
	int *Accepted;

} step_t;


static void worker_for_SS(void *_data, long i, int tid) // kt_for() callback
{
	step_t *step = (step_t*)_data;
	char *s = step->lines[i];
	int ReadLength =step->ReadLength;
	int DebugMode =step->DebugMode;
	int KmerSize =step->KmerSize;
	int IterationNo =step->IterationNo;
	int EditThreshold =step->EditThreshold;

	char RefSeq[ReadLength];
	char ReadSeq[ReadLength];

	int Checkpoint = -1;
    int rowNumber= -1;
    char *cigarSSWithoutMatch;

	//printf("%s\n",step->lines[i]);

	strncpy(ReadSeq, s, ReadLength);
	strncpy(RefSeq, s+ReadLength+1, ReadLength);

    /*
	for (n = 0; n < ReadLength; n++) {
		printf("%c",ReadSeq[n]);
	}
	printf("\t");
	for (n = 0; n < ReadLength; n++) {
		printf("%c",RefSeq[n]);
	}
	printf("\n");
    */

	step->Accepted[i] = SneakySnake(ReadLength, RefSeq, ReadSeq, EditThreshold, KmerSize, DebugMode, IterationNo, &Checkpoint, &rowNumber, &cigarSSWithoutMatch);
	//printf("main.c: i: %d TID:%d Accepted: %d\n",i, tid, step->Accepted[i]);

}

static void worker_for_Edlib(void *_data, long i, int tid) // kt_for() callback
{
    step_t *step = (step_t*)_data;
	char *s = step->lines[i];
	int ReadLength = step->ReadLength;

	int EditThreshold =step->EditThreshold;

    char RefSeq[ReadLength];
	char ReadSeq[ReadLength];

	strncpy(ReadSeq, s, ReadLength);
	strncpy(RefSeq, s+ReadLength+1, ReadLength);

    EdlibAlignResult resultEdlib ;
    char* cigarExtended ;
    char* cigarStandard ;
    int editDistance = 0;


    // EDLIB_MODE options: global (NW), prefix (SHW), infix (HW)
    resultEdlib = edlibAlign(RefSeq, ReadLength, ReadSeq, ReadLength, edlibNewAlignConfig(EditThreshold,EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

    if (resultEdlib.status == EDLIB_STATUS_OK && resultEdlib.editDistance != -1) {
        cigarExtended = edlibAlignmentToCigar(resultEdlib.alignment, resultEdlib.alignmentLength, EDLIB_CIGAR_EXTENDED);
        cigarStandard = edlibAlignmentToCigar(resultEdlib.alignment, resultEdlib.alignmentLength, EDLIB_CIGAR_STANDARD);

        /*
        printf("alignmentLength: %d\n", resultEdlib.alignmentLength);
        printf("startLocations[0]: %d\n", resultEdlib.startLocations[0]);
        printf("endLocations[0]: %d\n", resultEdlib.endLocations[0]);
        printf("alphabetLength: %d\n", resultEdlib.alphabetLength);
        printf("Cigar Extended format: %s\n",cigarExtended);
        printf("Cigar Standard format: %s\n",cigarStandard);
        */
        editDistance = resultEdlib.editDistance;

        Edlib_Edit[i] = editDistance;
        //printf("%d\t E:%d\n",i, editDistance);

        free(cigarExtended);
        free(cigarStandard);
    }
    if( resultEdlib.editDistance == -1) { // cannot find alignment having edit distance smaller than threshold
        editDistance = resultEdlib.editDistance;
        Edlib_Edit[i] = editDistance;
        //printf("%d\t E:%d\n",i, editDistance);
    }
    edlibFreeAlignResult(resultEdlib);

}

static void worker_for_SS_Edlib(void *_data, long i, int tid) // kt_for() callback
{
	step_t *step = (step_t*)_data;
	char *s = step->lines[i];
	int ReadLength =step->ReadLength;
	int DebugMode =step->DebugMode;
	int KmerSize =step->KmerSize;
	int IterationNo =step->IterationNo;
	int EditThreshold =step->EditThreshold;

	char RefSeq[ReadLength];
	char ReadSeq[ReadLength];

	int Checkpoint = -1;
    int rowNumber= -1;
    char *cigarSSWithoutMatch_;
    int dummy = 0;

	strncpy(ReadSeq, s, ReadLength);
	strncpy(RefSeq, s+ReadLength+1, ReadLength);

    int Runtime = 1;

    if (!Runtime) {

        EdlibAlignResult resultEdlibOnly;
        int EdlibAreaSize = -1;
        int fixedEdlibAreaSize = 5; //might change
        int startIndexSS = 0;
        int noOfMatch = 0;
        int EditThresholdCounter = EditThreshold;
        char* cigarExtended_ ;
        char* cigarStandard_;
        char *cigarTotalExtended = NULL;
        char *cigarTotalStandard = NULL;

        while( startIndexSS < ReadLength ){

            char *cigarSSExtended = NULL;
            char *cigarSSStandard = NULL;

            fprintf(stderr, "[[M::%s]] \n", __func__);

            //1. SneakySnake
            fprintf(stderr, "[M::%s] startIndexSS: %d\n", __func__,startIndexSS);
            char RefSeqSS[ReadLength-startIndexSS];
            char ReadSeqSS[ReadLength-startIndexSS];
            for (int i=startIndexSS;i<ReadLength;i++) {
                RefSeqSS[i-startIndexSS]=RefSeq[i];
                ReadSeqSS[i-startIndexSS]=ReadSeq[i];
            }
            fprintf(stderr, "[M::%s] EditThresholdCounter: %d\n", __func__,EditThresholdCounter);

            step->Accepted[i] = SneakySnake(ReadLength-startIndexSS, RefSeqSS, ReadSeqSS, EditThresholdCounter, ReadLength-startIndexSS, DebugMode, 0, &Checkpoint, &rowNumber, &cigarSSWithoutMatch_);
            fprintf(stderr, "[M::%s] i: %d TID:%d Accepted: %d\n", __func__,i, tid, step->Accepted[i]);

            if(step->Accepted[i] ==1) {
                fprintf(stderr, "[M::%s] Checkpoint: %d\n", __func__,Checkpoint);

                if ( EditThresholdCounter<0 ){
                    goto DONE;
                }

                // Decide the value of EdlibAreaSize
                if ( Checkpoint <= ((fixedEdlibAreaSize*2)+1) ){
                    EdlibAreaSize = Checkpoint;
                } else if ( Checkpoint > ((fixedEdlibAreaSize*2)+1) ) {
                    EdlibAreaSize = fixedEdlibAreaSize;
                }

                if( Checkpoint + EdlibAreaSize + startIndexSS >= ReadLength) {
                   EdlibAreaSize = -1;  //move on with SneakySnake
                }

                fprintf(stderr, "[M::%s] EdlibAreaSize: %d\n", __func__, EdlibAreaSize);
                fprintf(stderr, "[M::%s] cigarSSWithoutMatch_: %s\n", __func__,cigarSSWithoutMatch_);
                fprintf(stderr, "[M::%s] rowNumber: %d\n", __func__,rowNumber);

                char * cigarSSWithoutMatchStandard = NULL;
                char * cigarSSWithoutMatch = NULL;

                if (EdlibAreaSize == -1){ //move on with only SneakySnake
                    noOfMatch = Checkpoint ;

                    if (cigarSSWithoutMatch_ != NULL) {

                        cigarSSWithoutMatch = (char *)malloc( (strlen(cigarSSWithoutMatch_) + 1) * sizeof(char));
                        memcpy (cigarSSWithoutMatch, cigarSSWithoutMatch_, strlen(cigarSSWithoutMatch_));
                        cigarSSWithoutMatch [strlen(cigarSSWithoutMatch_)] = '\0';

                        cigarSSWithoutMatchStandard = (char *)malloc( (strlen(cigarSSWithoutMatch) + 1) * sizeof(char));
                        memcpy (cigarSSWithoutMatchStandard, cigarSSWithoutMatch, strlen(cigarSSWithoutMatch));
                        cigarSSWithoutMatchStandard [strlen(cigarSSWithoutMatch)] = '\0';

                        if(*(cigarSSWithoutMatchStandard+(strlen(cigarSSWithoutMatch)-1))== 'X') {
                            *(cigarSSWithoutMatchStandard+(strlen(cigarSSWithoutMatch)-1))= 'M';
                        }
                    }
                    //cigarSSWithoutMatch stays NULL
                    //cigarSSWithoutMatchStandard stays NULL


                } else {  //move on with only Edlib
                    noOfMatch = Checkpoint - EdlibAreaSize ;
                    // cigarSSWithoutMatchStandard is NULL
                    // cigarSSWithoutMatch is NULL
                }

                fprintf(stderr, "[M::%s] noOfMatch: %d\n", __func__, noOfMatch);
                fprintf(stderr, "[M::%s] cigarSSWithoutMatchExtended: %s\n", __func__, cigarSSWithoutMatch);
                fprintf(stderr, "[M::%s] cigarSSWithoutMatchStandard: %s\n", __func__, cigarSSWithoutMatchStandard);

                if(noOfMatch == 0 ) {
                    /*
                      when:
                        -Checkpoint = EdlibAreaSize, only Edlib
                        -Checkpoint = 0, only SneakySnake
                    */

                    if (cigarSSWithoutMatch != NULL) {
                        cigarSSExtended = (char *)malloc( sizeof(char) * ( strlen(cigarSSWithoutMatch) + 1));
                        memcpy (cigarSSExtended, cigarSSWithoutMatch, strlen(cigarSSWithoutMatch) );
                        cigarSSExtended [strlen(cigarSSWithoutMatch)] = '\0';
                    }
                    //else cigarSSExtended remains NULL
                    fprintf(stderr, "[M::%s] cigarSSExtended: %s\n", __func__,cigarSSExtended);

                    if (cigarSSWithoutMatchStandard != NULL) {
                         cigarSSStandard = (char *)malloc( sizeof(char) * (strlen(cigarSSWithoutMatchStandard) + 1));
                         memcpy (cigarSSStandard, cigarSSWithoutMatchStandard, strlen(cigarSSWithoutMatchStandard) );
                         cigarSSStandard [strlen(cigarSSWithoutMatchStandard)] = '\0';
                    }
                    //else cigarSSStandard remains NULL
                    fprintf(stderr, "[M::%s] cigarSSStandard: %s\n", __func__,cigarSSStandard);

                } else {
                    int counter = noOfMatch;
                    int noOfDigits = 0;
                    while (counter!=0) {
                        counter = counter/10;
                        noOfDigits++;
                    }

                    char *cigarSSWithMatch =  (char *)malloc( sizeof(char) * (noOfDigits + 1) );
                    char * number =  (char *)malloc( sizeof(char) * noOfDigits );
                    sprintf(number, "%d", noOfMatch);
                    memcpy (cigarSSWithMatch, number, noOfDigits );
                    cigarSSWithMatch [noOfDigits] = '\0';
                    free(number);

                    cigarSSWithMatch = (char *) realloc(cigarSSWithMatch, ( sizeof(char) * ( strlen(cigarSSWithMatch)  + 1 ) ) );
                    memcpy (  (cigarSSWithMatch + strlen(cigarSSWithMatch)) , "M",  2);
                    cigarSSWithMatch [strlen(cigarSSWithMatch)] = '\0';

                    fprintf(stderr, "[M::%s] cigarSSWithMatch: %s\n", __func__,cigarSSWithMatch);

                    if (cigarSSWithoutMatch != NULL) { //I-D var
                        cigarSSExtended = (char *)malloc( sizeof(char) * ( strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSExtended, cigarSSWithMatch, strlen(cigarSSWithMatch) );
                        cigarSSExtended [strlen(cigarSSWithMatch)] = '\0';

                        dummy = strlen(cigarSSExtended) + strlen(cigarSSWithoutMatch);
                        cigarSSExtended = (char *)realloc(cigarSSExtended, sizeof(char) * (strlen(cigarSSExtended) + strlen(cigarSSWithoutMatch) + 1) );
                        memcpy (cigarSSExtended + strlen(cigarSSExtended) , cigarSSWithoutMatch, strlen(cigarSSWithoutMatch) );
                        cigarSSExtended [dummy] = '\0';

                    } else {
                        cigarSSExtended = (char *) malloc( sizeof(char) * ( strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSExtended, cigarSSWithMatch, strlen(cigarSSWithMatch)  );
                        cigarSSExtended [strlen(cigarSSWithMatch)] = '\0';
                    }
                    fprintf(stderr, "[M::%s] cigarSSExtended: %s\n", __func__,cigarSSExtended);


                    if (cigarSSWithoutMatchStandard != NULL) {
                        cigarSSStandard = (char *)malloc( sizeof(char) * (strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSStandard, cigarSSWithMatch, strlen(cigarSSWithMatch) );
                        cigarSSStandard [strlen(cigarSSWithMatch)] = '\0';

                        dummy = strlen(cigarSSStandard) + strlen(cigarSSWithoutMatchStandard) ;
                        cigarSSStandard = (char *)realloc(cigarSSStandard, sizeof(char) * (strlen(cigarSSStandard) + strlen(cigarSSWithoutMatchStandard) + 1) );
                        memcpy (cigarSSStandard + strlen(cigarSSStandard) , cigarSSWithoutMatchStandard, strlen(cigarSSWithoutMatchStandard) );
                        cigarSSStandard [dummy] = '\0';

                    } else {
                        cigarSSStandard = (char *) malloc( sizeof(char) * (strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSStandard, cigarSSWithMatch, strlen(cigarSSWithMatch)  );
                        cigarSSStandard [strlen(cigarSSWithMatch)] = '\0';
                    }
                    fprintf(stderr, "[M::%s] cigarSSStandard: %s\n", __func__,cigarSSStandard);

                    free(cigarSSWithMatch);
                }

                if (cigarSSWithoutMatchStandard != NULL) {
                    free(cigarSSWithoutMatchStandard);
                }
                if (cigarSSWithoutMatch != NULL) {
                    free(cigarSSWithoutMatch);
                }

                // 2. Edlib
                if ( EdlibAreaSize != -1){ //Edlib function call

                    startIndexSS = startIndexSS + Checkpoint + EdlibAreaSize + 1;
                    int EdlibSequenceLength = (2 * EdlibAreaSize) + 1;

                    char RefSeqEdlib[EdlibSequenceLength];
                    char ReadSeqEdlib[EdlibSequenceLength];
                    for (int i=Checkpoint-EdlibAreaSize;i<=Checkpoint+EdlibAreaSize;i++){
                        RefSeqEdlib[i-(Checkpoint-EdlibAreaSize)]=RefSeq[i];
                        ReadSeqEdlib[i-(Checkpoint-EdlibAreaSize)]=ReadSeq[i];
                    }

                    // EDLIB_MODE options: global (NW), prefix (SHW), infix (HW)
                    resultEdlibOnly = edlibAlign(RefSeqEdlib, EdlibSequenceLength, ReadSeqEdlib, EdlibSequenceLength, edlibNewAlignConfig(EditThresholdCounter,
                        EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                    fprintf(stderr, "[M::%s] resultEdlibOnly.editDistance: %d\n", __func__,resultEdlibOnly.editDistance);

                    if (resultEdlibOnly.status == EDLIB_STATUS_OK && EditThresholdCounter - resultEdlibOnly.editDistance >= 0 && resultEdlibOnly.editDistance != -1) {
                        cigarExtended_ = edlibAlignmentToCigar(resultEdlibOnly.alignment, resultEdlibOnly.alignmentLength, EDLIB_CIGAR_EXTENDED);
                        cigarStandard_ = edlibAlignmentToCigar(resultEdlibOnly.alignment, resultEdlibOnly.alignmentLength, EDLIB_CIGAR_STANDARD);
                        fprintf(stderr, "[M::%s] EdlibArea-Cigar Extended format: %s\n", __func__,cigarExtended_);
                        fprintf(stderr, "[M::%s] EdlibArea-Cigar Standard format: %s\n", __func__,cigarStandard_);

                        EditThresholdCounter = EditThresholdCounter - resultEdlibOnly.editDistance;
                        fprintf(stderr, "[M::%s] EditThresholdCounter when EDLIB_STATUS_OK: %d\n", __func__, EditThresholdCounter);

                        if (cigarSSExtended == NULL) {
                            cigarSSExtended = (char *) malloc( sizeof(char) * (strlen(cigarExtended_) + 1) );
                            memcpy (cigarSSExtended, cigarExtended_, strlen(cigarExtended_) );
                            cigarSSExtended [strlen(cigarExtended_)] = '\0';
                        } else {
                            dummy = strlen(cigarExtended_) + strlen(cigarSSExtended);
                            cigarSSExtended = (char *) realloc(cigarSSExtended, sizeof(char) * (strlen(cigarExtended_) + strlen(cigarSSExtended) + 1));
                            memcpy (cigarSSExtended + strlen(cigarSSExtended) , cigarExtended_, strlen(cigarExtended_) );
                            cigarSSExtended [dummy] = '\0';
                        }
                        fprintf(stderr, "[M::%s] cigarSSExtended: %s\n", __func__, cigarSSExtended);

                        if (cigarTotalExtended == NULL ) {
                            cigarTotalExtended = (char *) malloc( sizeof(char) * (strlen(cigarSSExtended)+ 1));
                            memcpy (cigarTotalExtended, cigarSSExtended, strlen(cigarSSExtended) );
                            cigarTotalExtended [strlen(cigarSSExtended)] = '\0';
                        } else {
                            dummy = strlen(cigarTotalExtended) + strlen(cigarSSExtended) ;
                            cigarTotalExtended = (char *) realloc(cigarTotalExtended, sizeof(char) * (strlen(cigarTotalExtended) + strlen(cigarSSExtended) + 1));
                            memcpy (cigarTotalExtended + strlen(cigarTotalExtended) , cigarSSExtended, strlen(cigarSSExtended) );
                            cigarTotalExtended [dummy] = '\0';
                        }
                        fprintf(stderr, "[M::%s] cigarTotalExtended: %s\n", __func__, cigarTotalExtended);


                        if (cigarSSStandard == NULL) {
                            cigarSSStandard = (char *) malloc( sizeof(char) * (strlen(cigarStandard_) + 1 ));
                            memcpy (cigarSSStandard, cigarStandard_, strlen(cigarStandard_) );
                            cigarSSStandard [strlen(cigarStandard_)] = '\0';

                        } else {
                            dummy = strlen(cigarStandard_) + strlen(cigarSSStandard) ;
                            cigarSSStandard = (char *) realloc(cigarSSStandard, sizeof(char) * (strlen(cigarStandard_) + strlen(cigarSSStandard) + 1));
                            memcpy (cigarSSStandard + strlen(cigarSSStandard) , cigarStandard_, strlen(cigarStandard_) );
                            cigarSSStandard [dummy] = '\0';
                        }
                        fprintf(stderr, "[M::%s] cigarSSStandard: %s\n", __func__, cigarSSStandard);

                        if (cigarTotalStandard == NULL ) {
                            cigarTotalStandard = (char *) malloc( sizeof(char) * (strlen(cigarSSStandard) + 1 ));
                            memcpy (cigarTotalStandard, cigarSSStandard, strlen(cigarSSStandard) );
                            cigarTotalStandard [strlen(cigarSSStandard)] = '\0';
                        } else {
                            dummy = strlen(cigarTotalStandard) + strlen(cigarSSStandard);
                            cigarTotalStandard = (char *) realloc(cigarTotalStandard, sizeof(char) * (strlen(cigarTotalStandard) + strlen(cigarSSStandard) + 1 ));
                            memcpy (cigarTotalStandard + strlen(cigarTotalStandard) , cigarSSStandard, strlen(cigarSSStandard) );
                            cigarTotalStandard[dummy] = '\0';
                        }
                        fprintf(stderr, "[M::%s] cigarTotalStandard: %s\n", __func__, cigarTotalStandard);

                        free (cigarSSExtended);
                        free (cigarSSStandard);
                        free(cigarExtended_);
                        free(cigarStandard_);
                        edlibFreeAlignResult(resultEdlibOnly);

                    } else {
                        edlibFreeAlignResult(resultEdlibOnly);
                        goto DONE;
                    }

                } else { //SS function call (without Edlib)

                    startIndexSS = startIndexSS + Checkpoint + 1;
                    if (rowNumber == 0) {
                        EditThresholdCounter = EditThresholdCounter - 1;
                    } else {
                        EditThresholdCounter = EditThresholdCounter - rowNumber;
                    }
                    fprintf(stderr, "[M::%s] EditThresholdCounter after SS: %d\n", __func__, EditThresholdCounter);

                    if ( EditThresholdCounter < 0 ){
                        goto DONE;
                    }

                    if (cigarSSExtended != NULL) {

                        if(cigarTotalExtended == NULL) {
                              cigarTotalExtended = (char *)malloc( sizeof(char) * (strlen(cigarSSExtended) + 1 ));
                              memcpy (cigarTotalExtended, cigarSSExtended, strlen(cigarSSExtended) );
                              cigarTotalExtended [strlen(cigarSSExtended)] = '\0';

                        }else {
                             dummy = strlen(cigarSSExtended) + strlen(cigarTotalExtended) ;
                             cigarTotalExtended = (char *) realloc(cigarTotalExtended, sizeof(char) * (strlen(cigarSSExtended) + strlen(cigarTotalExtended) + 1));
                             memcpy (cigarTotalExtended + strlen(cigarTotalExtended) , cigarSSExtended, strlen(cigarSSExtended) );
                             cigarTotalExtended [dummy] = '\0';
                        }

                    }
                    fprintf(stderr, "[M::%s] cigarTotalExtended: %s\n", __func__, cigarTotalExtended);

                    if (cigarSSStandard != NULL) {

                        if(cigarTotalStandard == NULL) {
                              cigarTotalStandard = (char *)malloc( sizeof(char) * ( strlen(cigarSSStandard) + 1 ) );
                              memcpy (cigarTotalStandard, cigarSSStandard, strlen(cigarSSStandard) );
                              cigarTotalStandard [strlen(cigarSSStandard)] = '\0';
                        }else {
                             dummy = strlen(cigarTotalStandard) + strlen(cigarSSStandard) ;
                             cigarTotalStandard = (char *) realloc(cigarTotalStandard, sizeof(char) * (strlen(cigarTotalStandard) + strlen(cigarSSStandard) + 1));
                             memcpy (cigarTotalStandard + strlen(cigarTotalStandard) , cigarSSStandard, strlen(cigarSSStandard) );
                             cigarTotalStandard [dummy] = '\0';

                        }

                    }
                    fprintf(stderr, "[M::%s] cigarTotalStandard: %s\n", __func__, cigarTotalStandard);

                    if (cigarSSExtended !=NULL) {
                        free (cigarSSExtended);
                    }
                    if (cigarSSStandard !=NULL) {
                        free (cigarSSStandard);
                    }


                }
            }else {
                goto DONE;
            }
        }

        if ( EditThresholdCounter >= 0 ){
                    goto DONE_;
        }

        DONE:
        fprintf(stderr, "[M::%s] Done: EditThreshold exceed\n", __func__);
        fprintf(stderr, "[M::%s] startIndexSS supposed to be: %d\n", __func__,startIndexSS);
        free(cigarTotalExtended);
        free(cigarTotalStandard);
        //    printf("\n**Exonerate**\n\n");
        //    int r = system ("cd exonerate-2.2.0 && exonerate --model affine:global --exhaustive yes query.fasta target.fasta  --showvulgar no --showalignment yes --joinfilter 32 --showcigar yes --ryo \"%C\n\" ");
        //    printf("main.c: r: %d\n\n",r);
        //    int k = system ("cd ..");

        return;

        DONE_:
        fprintf(stderr, "[M::%s] Done_: EDLIB+SS is successful.\n", __func__);
        fprintf(stderr, "[M::%s] startIndexSS: %d\n", __func__,startIndexSS);
        printf("main.c: cigarTotalStandard: %s\n",cigarTotalStandard );
        printf("main.c: cigarTotalExtended: %s\n",cigarTotalExtended );

        free(cigarTotalExtended);
        free(cigarTotalStandard);

    }
    else { // no output

        EdlibAlignResult resultEdlibOnly;
        int EdlibAreaSize = -1;
        int fixedEdlibAreaSize = 5; //might change
        int startIndexSS = 0;
        int noOfMatch = 0;
        int EditThresholdCounter = EditThreshold;
        int editDistance = 0;
        char* cigarExtended_ ;
        char* cigarStandard_;
        char *cigarTotalExtended = NULL;
        char *cigarTotalStandard = NULL;

        while( startIndexSS < ReadLength ){

            char *cigarSSExtended = NULL;
            char *cigarSSStandard = NULL;

            //1. SneakySnake
            char RefSeqSS[ReadLength-startIndexSS];
            char ReadSeqSS[ReadLength-startIndexSS];
            for (int i=startIndexSS;i<ReadLength;i++) {
                RefSeqSS[i-startIndexSS]=RefSeq[i];
                ReadSeqSS[i-startIndexSS]=ReadSeq[i];
            }
            step->Accepted[i] = SneakySnake(ReadLength-startIndexSS, RefSeqSS, ReadSeqSS, EditThresholdCounter, ReadLength-startIndexSS, DebugMode, 0, &Checkpoint, &rowNumber, &cigarSSWithoutMatch_);

            if(step->Accepted[i] ==1) {

                if ( EditThresholdCounter<0 ){
                    goto DoNE;
                }

                // Decide the value of EdlibAreaSize
                if ( Checkpoint <= ((fixedEdlibAreaSize*2)+1) ){
                    EdlibAreaSize = Checkpoint;
                } else if ( Checkpoint > ((fixedEdlibAreaSize*2)+1) ) {
                    EdlibAreaSize = fixedEdlibAreaSize;
                }

                if( Checkpoint + EdlibAreaSize + startIndexSS >= ReadLength) {
                   EdlibAreaSize = -1;  //move on with SneakySnake
                }

                char * cigarSSWithoutMatchStandard = NULL;
                char * cigarSSWithoutMatch = NULL;

                if (EdlibAreaSize == -1){ //move on with only SneakySnake
                    noOfMatch = Checkpoint ;

                    if (cigarSSWithoutMatch_ != NULL) {

                        cigarSSWithoutMatch = (char *)malloc( (strlen(cigarSSWithoutMatch_) + 1) * sizeof(char));
                        memcpy (cigarSSWithoutMatch, cigarSSWithoutMatch_, strlen(cigarSSWithoutMatch_));
                        cigarSSWithoutMatch [strlen(cigarSSWithoutMatch_)] = '\0';

                        cigarSSWithoutMatchStandard = (char *)malloc( (strlen(cigarSSWithoutMatch) + 1) * sizeof(char));
                        memcpy (cigarSSWithoutMatchStandard, cigarSSWithoutMatch, strlen(cigarSSWithoutMatch));
                        cigarSSWithoutMatchStandard [strlen(cigarSSWithoutMatch)] = '\0';

                        if(*(cigarSSWithoutMatchStandard+(strlen(cigarSSWithoutMatch)-1))== 'X') {
                            *(cigarSSWithoutMatchStandard+(strlen(cigarSSWithoutMatch)-1))= 'M';
                        }
                    }
                    //cigarSSWithoutMatch stays NULL
                    //cigarSSWithoutMatchStandard stays NULL


                } else {  //move on with only Edlib
                    noOfMatch = Checkpoint - EdlibAreaSize ;
                    // cigarSSWithoutMatchStandard is NULL
                    // cigarSSWithoutMatch is NULL
                }

                if(noOfMatch == 0 ) {
                    /*
                      when:
                        -Checkpoint = EdlibAreaSize, only Edlib
                        -Checkpoint = 0, only SneakySnake
                    */

                    if (cigarSSWithoutMatch != NULL) {
                        cigarSSExtended = (char *)malloc( sizeof(char) * ( strlen(cigarSSWithoutMatch) + 1));
                        memcpy (cigarSSExtended, cigarSSWithoutMatch, strlen(cigarSSWithoutMatch) );
                        cigarSSExtended [strlen(cigarSSWithoutMatch)] = '\0';
                    }
                    //else cigarSSExtended remains NULL

                    if (cigarSSWithoutMatchStandard != NULL) {
                         cigarSSStandard = (char *)malloc( sizeof(char) * (strlen(cigarSSWithoutMatchStandard) + 1));
                         memcpy (cigarSSStandard, cigarSSWithoutMatchStandard, strlen(cigarSSWithoutMatchStandard) );
                         cigarSSStandard [strlen(cigarSSWithoutMatchStandard)] = '\0';
                    }
                    //else cigarSSStandard remains NULL

                } else {
                    int counter = noOfMatch;
                    int noOfDigits = 0;
                    while (counter!=0) {
                        counter = counter/10;
                        noOfDigits++;
                    }

                    char *cigarSSWithMatch =  (char *)malloc( sizeof(char) * (noOfDigits + 1) );
                    char * number =  (char *)malloc( sizeof(char) * noOfDigits );
                    sprintf(number, "%d", noOfMatch);
                    memcpy (cigarSSWithMatch, number, noOfDigits );
                    cigarSSWithMatch [noOfDigits] = '\0';
                    free(number);

                    cigarSSWithMatch = (char *) realloc(cigarSSWithMatch, ( sizeof(char) * ( strlen(cigarSSWithMatch)  + 1 ) ) );
                    memcpy (  (cigarSSWithMatch + strlen(cigarSSWithMatch)) , "M",  2);
                    cigarSSWithMatch [strlen(cigarSSWithMatch)] = '\0';

                    if (cigarSSWithoutMatch != NULL) { //I-D var
                        cigarSSExtended = (char *)malloc( sizeof(char) * ( strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSExtended, cigarSSWithMatch, strlen(cigarSSWithMatch) );
                        cigarSSExtended [strlen(cigarSSWithMatch)] = '\0';

                        dummy = strlen(cigarSSExtended) + strlen(cigarSSWithoutMatch);
                        cigarSSExtended = (char *)realloc(cigarSSExtended, sizeof(char) * (strlen(cigarSSExtended) + strlen(cigarSSWithoutMatch) + 1) );
                        memcpy (cigarSSExtended + strlen(cigarSSExtended) , cigarSSWithoutMatch, strlen(cigarSSWithoutMatch) );
                        cigarSSExtended [dummy] = '\0';

                    } else {
                        cigarSSExtended = (char *) malloc( sizeof(char) * ( strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSExtended, cigarSSWithMatch, strlen(cigarSSWithMatch)  );
                        cigarSSExtended [strlen(cigarSSWithMatch)] = '\0';
                    }

                    if (cigarSSWithoutMatchStandard != NULL) {
                        cigarSSStandard = (char *)malloc( sizeof(char) * (strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSStandard, cigarSSWithMatch, strlen(cigarSSWithMatch) );
                        cigarSSStandard [strlen(cigarSSWithMatch)] = '\0';

                        dummy = strlen(cigarSSStandard) + strlen(cigarSSWithoutMatchStandard) ;
                        cigarSSStandard = (char *)realloc(cigarSSStandard, sizeof(char) * (strlen(cigarSSStandard) + strlen(cigarSSWithoutMatchStandard) + 1) );
                        memcpy (cigarSSStandard + strlen(cigarSSStandard) , cigarSSWithoutMatchStandard, strlen(cigarSSWithoutMatchStandard) );
                        cigarSSStandard [dummy] = '\0';

                    } else {
                        cigarSSStandard = (char *) malloc( sizeof(char) * (strlen(cigarSSWithMatch) + 1) );
                        memcpy (cigarSSStandard, cigarSSWithMatch, strlen(cigarSSWithMatch)  );
                        cigarSSStandard [strlen(cigarSSWithMatch)] = '\0';
                    }

                    free(cigarSSWithMatch);
                }

                if (cigarSSWithoutMatchStandard != NULL) {
                    free(cigarSSWithoutMatchStandard);
                }
                if (cigarSSWithoutMatch != NULL) {
                    free(cigarSSWithoutMatch);
                }

                // 2. Edlib
                if ( EdlibAreaSize != -1){ //Edlib function call

                    startIndexSS = startIndexSS + Checkpoint + EdlibAreaSize + 1;
                    int EdlibSequenceLength = (2 * EdlibAreaSize) + 1;

                    char RefSeqEdlib[EdlibSequenceLength];
                    char ReadSeqEdlib[EdlibSequenceLength];
                    for (int i=Checkpoint-EdlibAreaSize;i<=Checkpoint+EdlibAreaSize;i++){
                        RefSeqEdlib[i-(Checkpoint-EdlibAreaSize)]=RefSeq[i];
                        ReadSeqEdlib[i-(Checkpoint-EdlibAreaSize)]=ReadSeq[i];
                    }

                    // EDLIB_MODE options: global (NW), prefix (SHW), infix (HW)
                    resultEdlibOnly = edlibAlign(RefSeqEdlib, EdlibSequenceLength, ReadSeqEdlib, EdlibSequenceLength, edlibNewAlignConfig(EditThresholdCounter,
                        EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

                    if (resultEdlibOnly.status == EDLIB_STATUS_OK && EditThresholdCounter - resultEdlibOnly.editDistance >= 0 && resultEdlibOnly.editDistance != -1) {
                        cigarExtended_ = edlibAlignmentToCigar(resultEdlibOnly.alignment, resultEdlibOnly.alignmentLength, EDLIB_CIGAR_EXTENDED);
                        cigarStandard_ = edlibAlignmentToCigar(resultEdlibOnly.alignment, resultEdlibOnly.alignmentLength, EDLIB_CIGAR_STANDARD);

                        EditThresholdCounter = EditThresholdCounter - resultEdlibOnly.editDistance;

                        if (cigarSSExtended == NULL) {
                            cigarSSExtended = (char *) malloc( sizeof(char) * (strlen(cigarExtended_) + 1) );
                            memcpy (cigarSSExtended, cigarExtended_, strlen(cigarExtended_) );
                            cigarSSExtended [strlen(cigarExtended_)] = '\0';
                        } else {
                            dummy = strlen(cigarExtended_) + strlen(cigarSSExtended);
                            cigarSSExtended = (char *) realloc(cigarSSExtended, sizeof(char) * (strlen(cigarExtended_) + strlen(cigarSSExtended) + 1));
                            memcpy (cigarSSExtended + strlen(cigarSSExtended) , cigarExtended_, strlen(cigarExtended_) );
                            cigarSSExtended [dummy] = '\0';
                        }

                        if (cigarTotalExtended == NULL ) {
                            cigarTotalExtended = (char *) malloc( sizeof(char) * (strlen(cigarSSExtended)+ 1));
                            memcpy (cigarTotalExtended, cigarSSExtended, strlen(cigarSSExtended) );
                            cigarTotalExtended [strlen(cigarSSExtended)] = '\0';
                        } else {
                            dummy = strlen(cigarTotalExtended) + strlen(cigarSSExtended) ;
                            cigarTotalExtended = (char *) realloc(cigarTotalExtended, sizeof(char) * (strlen(cigarTotalExtended) + strlen(cigarSSExtended) + 1));
                            memcpy (cigarTotalExtended + strlen(cigarTotalExtended) , cigarSSExtended, strlen(cigarSSExtended) );
                            cigarTotalExtended [dummy] = '\0';
                        }

                        if (cigarSSStandard == NULL) {
                            cigarSSStandard = (char *) malloc( sizeof(char) * (strlen(cigarStandard_) + 1 ));
                            memcpy (cigarSSStandard, cigarStandard_, strlen(cigarStandard_) );
                            cigarSSStandard [strlen(cigarStandard_)] = '\0';

                        } else {
                            dummy = strlen(cigarStandard_) + strlen(cigarSSStandard) ;
                            cigarSSStandard = (char *) realloc(cigarSSStandard, sizeof(char) * (strlen(cigarStandard_) + strlen(cigarSSStandard) + 1));
                            memcpy (cigarSSStandard + strlen(cigarSSStandard) , cigarStandard_, strlen(cigarStandard_) );
                            cigarSSStandard [dummy] = '\0';
                        }

                        if (cigarTotalStandard == NULL ) {
                            cigarTotalStandard = (char *) malloc( sizeof(char) * (strlen(cigarSSStandard) + 1 ));
                            memcpy (cigarTotalStandard, cigarSSStandard, strlen(cigarSSStandard) );
                            cigarTotalStandard [strlen(cigarSSStandard)] = '\0';
                        } else {
                            dummy = strlen(cigarTotalStandard) + strlen(cigarSSStandard);
                            cigarTotalStandard = (char *) realloc(cigarTotalStandard, sizeof(char) * (strlen(cigarTotalStandard) + strlen(cigarSSStandard) + 1 ));
                            memcpy (cigarTotalStandard + strlen(cigarTotalStandard) , cigarSSStandard, strlen(cigarSSStandard) );
                            cigarTotalStandard[dummy] = '\0';
                        }

                        free (cigarSSExtended);
                        free (cigarSSStandard);
                        free(cigarExtended_);
                        free(cigarStandard_);
                        edlibFreeAlignResult(resultEdlibOnly);

                    } else {
                        edlibFreeAlignResult(resultEdlibOnly);
                        goto DoNE;
                    }

                } else { //SS function call (without Edlib)

                    startIndexSS = startIndexSS + Checkpoint + 1;
                    if (rowNumber == 0) {
                        EditThresholdCounter = EditThresholdCounter - 1;
                    } else {
                        EditThresholdCounter = EditThresholdCounter - rowNumber;
                    }

                    if ( EditThresholdCounter < 0 ){
                        goto DoNE;
                    }

                    if (cigarSSExtended != NULL) {

                        if(cigarTotalExtended == NULL) {
                              cigarTotalExtended = (char *)malloc( sizeof(char) * (strlen(cigarSSExtended) + 1 ));
                              memcpy (cigarTotalExtended, cigarSSExtended, strlen(cigarSSExtended) );
                              cigarTotalExtended [strlen(cigarSSExtended)] = '\0';

                        }else {
                             dummy = strlen(cigarSSExtended) + strlen(cigarTotalExtended) ;
                             cigarTotalExtended = (char *) realloc(cigarTotalExtended, sizeof(char) * (strlen(cigarSSExtended) + strlen(cigarTotalExtended) + 1));
                             memcpy (cigarTotalExtended + strlen(cigarTotalExtended) , cigarSSExtended, strlen(cigarSSExtended) );
                             cigarTotalExtended [dummy] = '\0';
                        }

                    }

                    if (cigarSSStandard != NULL) {

                        if(cigarTotalStandard == NULL) {
                              cigarTotalStandard = (char *)malloc( sizeof(char) * ( strlen(cigarSSStandard) + 1 ) );
                              memcpy (cigarTotalStandard, cigarSSStandard, strlen(cigarSSStandard) );
                              cigarTotalStandard [strlen(cigarSSStandard)] = '\0';
                        }else {
                             dummy = strlen(cigarTotalStandard) + strlen(cigarSSStandard) ;
                             cigarTotalStandard = (char *) realloc(cigarTotalStandard, sizeof(char) * (strlen(cigarTotalStandard) + strlen(cigarSSStandard) + 1));
                             memcpy (cigarTotalStandard + strlen(cigarTotalStandard) , cigarSSStandard, strlen(cigarSSStandard) );
                             cigarTotalStandard [dummy] = '\0';

                        }

                    }

                    if (cigarSSExtended !=NULL) {
                        free (cigarSSExtended);
                    }
                    if (cigarSSStandard !=NULL) {
                        free (cigarSSStandard);
                    }

                }
            }else {
                goto DoNE;
            }
        }

        if ( EditThresholdCounter >= 0 ){
                    goto DoNE_;
        }

        DoNE:
        editDistance = -1 ;
        //printf("%d\t ES:%d\n",i, editDistance);
        SS_Edlib_Edit[i] = editDistance;

        free(cigarTotalExtended);
        free(cigarTotalStandard);
        //    printf("\n**Exonerate**\n\n");
        //    int r = system ("cd exonerate-2.2.0 && exonerate --model affine:global --exhaustive yes query.fasta target.fasta  --showvulgar no --showalignment yes --joinfilter 32 --showcigar yes --ryo \"%C\n\" ");
        //    printf("main.c: r: %d\n\n",r);
        //    int k = system ("cd ..");

        return;

        DoNE_:
        editDistance = EditThreshold - EditThresholdCounter ;
        //printf("%d\t ES:%d\n",i, editDistance);
        SS_Edlib_Edit[i] = editDistance;

        //printf("main.c: cigarTotalStandard: %s\n",cigarTotalStandard );
        //printf("main.c: cigarTotalExtended: %s\n",cigarTotalExtended );

        free(cigarTotalExtended);
        free(cigarTotalStandard);

    }


}

int main(int argc, const char * const argv[]) {
	
	if (argc!=9){
		fprintf(stderr, "missing argument..\n./main [DebugMode] [KmerSize] [ReadLength] [IterationNo] [ReadRefFile] [# of reads] [# of threads] [EditThreshold]\n");
		exit(-1);
	}

	FILE * fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	int i;
	int loopPar;
	long start, end;
    struct timeval timecheck;
	int threads = atoi(argv[7]);
	
	printf("This system has %d processors configured and %d processors available.\n", get_nprocs_conf(), get_nprocs());
	
	fp = fopen(argv[5], "r");
	if (!fp){
		fprintf(stderr, "Sorry, the file does not exist or you do not have access permission\n");
	}
	else {
		
        if (threads==1) {
            int ReadLength = atoi(argv[3]);
            int DebugMode = atoi(argv[1]);
            int KmerSize = atoi(argv[2]);
            int IterationNo = atoi(argv[4]);
            int EditThreshold = atoi(argv[8]);
            int Accepted =0;
            char RefSeq[ReadLength];
            char ReadSeq[ReadLength];

            int Checkpoint = -1;
            int rowNumber = -1;
            char *cigarSSWithoutMatch;

            for (i = 0; i < atoi(argv[6]); i++) {
                read = getline(&line, &len, fp);
                char *s = strdup(line);
                strncpy(ReadSeq, s, ReadLength);
                strncpy(RefSeq, s+ReadLength+1, ReadLength);

                /*
                for (n = 0; n < ReadLength; n++) {
                printf("%c",ReadSeq[n]);
                }
                printf("\t");
                for (n = 0; n < ReadLength; n++) {
                printf("%c",RefSeq[n]);
                }
                printf("\n");
                */

				if (SneakySnake(ReadLength, RefSeq, ReadSeq, EditThreshold, KmerSize, DebugMode, IterationNo,&Checkpoint, &rowNumber, &cigarSSWithoutMatch ))
					Accepted++;
				//printf("i: %d TID:%d Accepted: %d\n",i, tid, step->Accepted[i]);
                    
            }
             printf("Data: %s\tThreads: %d\tE: %d\tAccepted: %d\tRejected:%d\n", argv[5], threads, EditThreshold, Accepted, atoi(argv[6])-Accepted);
        } else { //when threads >1
            step_t *s;
            s = calloc(1, sizeof(step_t));
            s->lines = calloc(atoi(argv[6]), sizeof(char*));
            s->Accepted = calloc(atoi(argv[6]), sizeof(int));
            s->DebugMode = atoi(argv[1]);
            s->KmerSize = atoi(argv[2]);
            s->ReadLength  = atoi(argv[3]);
            s->IterationNo  = atoi(argv[4]);
            s->EditThreshold = atoi(argv[8]);

            Edlib_Edit = (int *)malloc( sizeof(int) * (atoi(argv[6])));
            SS_Edlib_Edit = (int *)malloc( sizeof(int) * (atoi(argv[6])));

            int difference = 0;
            int positive = 0;
            int neutral = 0;
            int negative = 0;

            for (i = 0; i < atoi(argv[6]); i++) {
                read = getline(&line, &len, fp);
                //printf("%d",s->n_lines);
                s->lines[s->n_lines] = strdup(line);
                s->Accepted[s->n_lines] = 0;
                //free(s->lines[s->n_lines]);
                ++s->n_lines;	
            }

            // Comparison
            gettimeofday(&timecheck, NULL);
            start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

            for (loopPar = 3; loopPar<=3;loopPar++) {
                kt_for(threads, worker_for_Edlib, s, atoi(argv[6]) );
                kt_for(threads, worker_for_SS_Edlib, s, atoi(argv[6]));
            }

            gettimeofday(&timecheck, NULL);
            end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

            printf("\tSneaky Snake and Edlib Total Time: %ld milliseconds\n", (end - start));

            // Save the outputs to a file
            FILE *f = fopen("File.txt", "w");
            if (f == NULL)  {
                printf("Error!\n");
                exit(1);
            }
            printf("%6s  %6s %6s %6s\n", "Seq No","X" , "Y" , "X-Y" );
            fprintf(f,"%6s  %6s %6s %6s\n", "Seq No","X" , "Y" , "X-Y" );

            for (int x = 0; x < atoi(argv[6]); x++) {
                if ( Edlib_Edit[x] == -1 &&  SS_Edlib_Edit[x] != -1) {
                    difference = SS_Edlib_Edit[x] -  Edlib_Edit[x] ;
                }else if ( Edlib_Edit[x] != -1 &&  SS_Edlib_Edit[x] == -1){
                    difference = SS_Edlib_Edit[x] -  Edlib_Edit[x] ;
                } else if ( Edlib_Edit[x] == -1 &&  SS_Edlib_Edit[x] == -1) {
                    difference= -1;
                }else{
                    difference = Edlib_Edit[x] - SS_Edlib_Edit[x]  ;
                }
                if ( difference > 0) {
                    positive++;
                } else if (difference == 0) {
                    neutral++;
                }
               //printf("%6d %6d %6d %6d\n", x, Edlib_Edit[x], SS_Edlib_Edit[x], difference );
                fprintf(f,"%6d %6d %6d %6d\n", x, Edlib_Edit[x], SS_Edlib_Edit[x], difference );
            }
            negative = atoi(argv[6]) - neutral - positive;
            printf("For %s\n #positive %d\n #neutral %d\n #negative %d\n",argv[5] ,positive,neutral,negative);
            fprintf(f,"For %s\n #positive %d\n #neutral %d\n #negative %d\n",argv[5] ,positive,neutral,negative);

            fclose(f);


            for (i = 0; i < atoi(argv[6]); i++) {
                free(s->lines[i]);
                ++s->n_lines;
            } 
            free(s->lines); free(s->Accepted); free(s);
        }
	}

	fclose(fp);
	return 0;
}
