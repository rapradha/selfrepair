<?php

/* first step is to determine if the input graph does not have any leaves.
 For this, check the adjacency matrix's rows, the sum of each row should be >= 2.
 */

 function AnyLeaves($AdjMatrix, $n){
 	$leafexists = false;

 	for($i=0;$i<$n;$i++){
 		$rowSum = 0;
 		for($j=0;$j<$n;$j++){
 			$rowSum += $AdjMatrix[$i][$j];
 		}//end for $j<$nc

 		if($rowSum < 2){
 			$leafexists = true;
 			return $leafexists;
 		}
 	}//end for $i<$n

 	return $leafexists;
 }//end function AnyLeaves

 function ReadAdjMatrix(){
 	$IPArray = array();
 	$ipfilename = "AdjMatrix.txt";
 	$IPArray = split("\n", file_get_contents($ipfilename));

 	return $IPArray;
 }//end function ReadAdjMatrix

 function GetInMatrixSize($AdjMatrix){
 	return count($AdjMatrix);
 }//end function GetInMatrixSize

 function AllocateSpace($mat, $n){
 	/* allocate space */
 	for($i=0;$i<$n;$i++){
 		for($j=0;$j<$n;$j++){
 			$mat[$i][$j] = 0;
 		}
 	}
 	return $mat;
 }//end allocateSpace

 function convertInMatrixToAdjMatrix($InMatrix, $n){
 	$convMatrix = array();
 	$convMatrix = AllocateSpace($convMatrix,$n);
 	
 	for($i=0;$i<$n;$i++){
 		$str = $InMatrix[$i];
 		$arVal = array();

 		$arVal = split(" ", $str);

 		for($j=0;$j<$n;$j++){
 			$convMatrix[$i][$j] = (int)$arVal[$j];
 		}
 	}
 	return $convMatrix;
 }//end convertInMatrixToAdjMatrix

function getOpenNeighborhoodInfo($AdjMatrix,$n){
	$OpenNeighbors = array();
	
	for($i=0;$i<$n;$i++){
		$OpenNeighbors[$i] = array();
	}

	for($i=0;$i<$n;$i++){
		for($j=0;$j<$n;$j++){
			if($AdjMatrix[$i][$j] ==  1){
				array_push($OpenNeighbors[$i], $j);
			}
		}
	}

	return $OpenNeighbors;
}//end getOpenNeighborhoodInfo

/****************************************************************************************************/
function getNextPermutation( $aPermutableItems, $iPermutationSize, $aPreviousPermutation = NULL )
{
  $aNextPermutation = $aPreviousPermutation;
  $iLastIndex       = $iPermutationSize - 1;
  $iPermutableItems = count($aPermutableItems);

  // Any previous permutation ?
  if( $aPreviousPermutation )
  {
    // Loop the elements backwards
    for( $i = $iLastIndex; $i >= 0; $i-- )
    {
      // Can the current element be incremented without reaching the limit ?
      if( ++$aNextPermutation[ $i ] >= $iPermutableItems )
      {
        // Increment the previous element
        $iPrevValue = ++$aNextPermutation[ $i - 1 ];
        // Reset the current element with the value of the previous plus one
        $iNextValue = $aNextPermutation[ $i ] = $iPrevValue + 1;
        // Skip the previous element because it was just incremented
        $i--;
        // If one of the two elements reached the limit, we are in the exit condition
        if( $iPrevValue >= $iPermutableItems || $iNextValue >= $iPermutableItems )
          return FALSE;
      }
      // Limit still to be reached for the i-th element, skip previous ones
      else
        break;
    }
  }
  // Am i able to generate the first permutation ?
  else if( $iPermutationSize <= $iPermutableItems )
    $aNextPermutation = range( 0, $iLastIndex );
  // Permutation impossible to generate because we don't have enough elements in the main set
  else
    $aNextPermutation = FALSE;

  return $aNextPermutation;
}

function computeSize($cn){
  $numerator = 1;
  $denominator1 = 1;
  $denominator2 = 2;

  for($i=$cn;$i>1;$i--){
    $numerator *= $i;
  }

  for($j=($cn-2);$j>1;$j--){
    $denominator1 *= $j;
  }

  return ($numerator/($denominator1*$denominator2));
}

function getAllTwoCombinationOfOpenNeighbors($aSize){
  $iPerm  = 0;
  $aPrev  = NULL;
  $aItems = array();

  for($i=0;$i<$aSize;$i++){
    $aItems[$i] = $i;
  }

  $permArray = array();
  $permArraySize = computeSize($aSize);
  for($i=0;$i<$permArraySize;$i++){
    $permArray[$i] = array();
  }

  $count = 0;
  while( ($aPrev = getNextPermutation( $aItems, 2, $aPrev )) != FALSE )
  {
    $temp = array();
    $temp = $aPrev;
    
    array_push($permArray[$count], $temp[0]);
    array_push($permArray[$count], $temp[1]);
    $count++;
  }

  // print_r($permArray);

  return $permArray;

}

/****************************************************************************************************/

function checkForThreeCycle($AdjMatrix,$Neighbors, $n){
	/* first check for three cycle */
	$ThreeCycleVertices = array(); //this keeps track of vertices that together with all its neighbors forms three cycles
	$FourCycleVertices = array();

	/* take one vertex and its neighboring vertices in pairs,
	if these three vertices share an edge between each other, then there is a three cycle */

	for($vertex=0;$vertex<$n;$vertex++){
		/* i has the first vertex */

		/* check if all neigbors form a complete graph */
		$edgeList = array();
		$numNeighbors = count($Neighbors[$vertex]);

		// print("Start ---\n Showing neighbors first\n");
		// print_r($Neighbors[$vertex]);
		// print("vertex $vertex no. of neighbors: $numNeighbors\n");
		$edgeListSize = computeSize($numNeighbors);
		// print("edge list size: $edgeListSize\n");
		$TotalNumOfEdges = $numNeighbors*($numNeighbors - 1)/2;
		// print("total num edges: $TotalNumOfEdges\n");
		
		if($TotalNumOfEdges == 1){
			$n1 = $Neighbors[$vertex][0];
	  	$n2 = $Neighbors[$vertex][1];

	  	// print("one edge neighbors: $nbr1 $nbr2\n");
	 		if($AdjMatrix[$n1][$n2] == 1){
	 			if(!in_array($vertex, $ThreeCycleVertices)){
	  				array_push($ThreeCycleVertices, $vertex);	
	  		}
	  	}else{
	  			/* do these neighbors belong to a four cycle? */
	  			$nbr1ON = array();
	  			$nbr2ON = array();
	  			$commonNeigbors = array();

	  			$nbr1ON = $Neighbors[$n1];
	  			$nbr2ON = $Neighbors[$n2];

	  			/* check if any of nbr1's  neighbor is also nbr2's neighbor, apart from vertex */
	  			for($j=0; $j < count($nbr1ON); $j++){
	  				$nbr = $nbr1ON[$j];

	  				if(in_array($nbr, $nbr2ON) && $nbr != $vertex){
	  					array_push($commonNeigbors, $nbr);
	  				}
	  			}

	  			if(count($commonNeigbors) > 0){
	  				/* the neighbors belong to a four cycle */
	  				array_push($FourCycleVertices, $vertex);
	  			}
	  		}//end if-else 3-4 cycle
		}else{
			for($i=0;$i<$edgeListSize;$i++){
	    	$edgeList[$i] = array();
	  	}

	  	$edgeList = getAllTwoCombinationOfOpenNeighbors($edgeListSize);
	  	// print("edge list\n");
	  	// print_r($edgeList);

	  	/* now for each of these two combinations, check if the vertices share an edge */
	  	for($i=0;$i<$edgeListSize;$i++){
	  		$nbr1 = $edgeList[$i][0];
	  		$nbr2 = $edgeList[$i][1];

	  		$n1 = $Neighbors[$vertex][$nbr1];
	  		$n2 = $Neighbors[$vertex][$nbr2];

	  		// print("neighbors: $nbr1 $nbr2\n");
	  		if($AdjMatrix[$n1][$n2] == 1){
	  			if(!in_array($vertex, $ThreeCycleVertices)){
	  				array_push($ThreeCycleVertices, $vertex);	
	  			}
	  		}else{
	  			/* do these neighbors belong to a four cycle? */
	  			$nbr1ON = array();
	  			$nbr2ON = array();
	  			$commonFourNeigbors = array();

	  			$nbr1ON = $Neighbors[$n1];
	  			$nbr2ON = $Neighbors[$n2];

	  			/* check if any of nbr1's  neighbor is also nbr2's neighbor, apart from vertex */
	  			for($j=0; $j < count($nbr1ON); $j++){
	  				$nbr = $nbr1ON[$j];

	  				if(in_array($nbr, $nbr2ON) && $nbr != $vertex){
	  					array_push($commonFourNeigbors, $nbr);
	  				}
	  			}

	  			if(count($commonFourNeigbors) > 0){
	  				/* the neighbors belong to a four cycle */
	  				if(!in_array($vertex, $FourCycleVertices)){
	  					array_push($FourCycleVertices, $vertex);	
	  				}
	  			}
	  		}//end if-else 3-4 cycle
	  	}//end for i<edgeListSize 
		}//end if-else TotalNumOfEdges == 1
	}//end for vertex<n

	print("all three cycle vertices\n");
	print_r($ThreeCycleVertices);

	print("all four cycle vertices\n");
	print_r($FourCycleVertices);

	/* combine both the three cycles and fourcycles, remove duplicates, and check if all vertices are present.  If yes, then graph is self-repairing else not self-repairing. */
	$cycleVertices = array();

	$cycleVertices = array_unique( array_merge($ThreeCycleVertices, $FourCycleVertices));
	$checkSum = ($n-1)*($n/2);
	$aSum = array_sum($cycleVertices);

	if($checkSum == $aSum){
		return true;
	}else{
		return false;
	}
}//end checkForThreeCycle

function SelfRepair(){
 	$InMatrix = array();
 	$AdjMatrix = array();
 	$Neighbors = array();
 	$ThreeCycleVertex = array();
 	$FourCycleVertex = array();

 	$n = 0;

 	$InMatrix = ReadAdjMatrix();
 	$n = GetInMatrixSize($InMatrix);
 	$AdjMatrix = convertInMatrixToAdjMatrix($InMatrix, $n);

 	/* Step 1: check if there are any leaves */
 	if(AnyLeaves($AdjMatrix,$n)){
 		print("has leaves, cannot self-repair\n");
 		return;
 	}else{
 		print("no leaves, check further for self-repair\n");

 		for($i=0;$i<$n;$i++){
 			$Neighbors[$i] = array();
 		}

 		/* Step 2: For each vertex gather info about its open neighbors, neighbors at one step away */
 		$Neighbors = getOpenNeighborhoodInfo($AdjMatrix,$n);

 		/* Step 3: For each vertex v and for all neighbors in its open neighborhood pair-wise, check if the three vertices form a three cycle or if they belong to a four cycle */
 		if(checkForThreeCycle($AdjMatrix,$Neighbors,$n)){
 			print("self repair possible\n");
 		}else{
 			print("self repair not possible\n");
 		}
 	}//end if-else no leaves

 	return;
 }//end function self-repair
 
 SelfRepair();
?>