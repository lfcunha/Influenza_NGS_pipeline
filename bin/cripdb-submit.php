#!/usr/bin/env php

<?php


$options = array(
                    PDO::MYSQL_ATTR_INIT_COMMAND => 'SET NAMES utf8',
                    PDO::FETCH_ASSOC
                );


/*  read db configuration and set up connection */

$db_config="/hpc/users/cunhal01/my.db.php.conf";
if(!file_exists($db_config)) die("Please provide db configuration file\n");
$db_param = unserialize(fread(fopen($db_config, "r"), filesize($db_config)));
$dbh= new \PDO('mysql:host='.$db_param["host"].';dbname='.$db_param["db"], $db_param['user'], $db_param['passw'], $options);


if ($argv[1]!="0") {
  
    //** get extractID of the sample (from joined Isolates + Extracts view) **//
    $query="SELECT idExtract FROM `ExtractsIsolates` WHERE `Sample_Identifier` LIKE '" . $argv[1] ."'" ;
    $sth=$dbh->prepare($query);
    $sth->execute();
    $result=$sth->fetch(PDO::FETCH_ASSOC);
    $extractID=intval($result["idExtract"]);
    //****//



    //** read the fasta sequence and build a hash of the segment's sequences **// 
    $filename=$argv[1]."_curated.fa";
    if(!file_exists($filename)) die("Extract does not exits - check name\n");


    /*
    * Process Sequences table
    */

    /*
    * read and encode intact fasta file to upload to 'assemblies' table
    */
    $fasta = base64_encode(fread(fopen($filename, "r"), filesize($filename)));
    /*
    * open fasta file to process each segment to upload to 'Sequences' table
    */
	$handle = fopen($filename, "r");
	$complete=1;
	if ($handle) {
        $count=array();
        $segments=array();
        $sequence="";
        $line = fgets($handle);
        if (strpos($line, ">")>-1){
            if(!strpos($line, "complete")){
                $complete=0;
                echo $complete;
            };
	    if (array_key_exists(intval($line[1]), $count)){
		$count[intval($line[1])]+=1;
	    }
	    else {
		$count[intval($line[1])]=1;
	    }

            $segment=$line[1];
            $info=explode("|", $line);
	    $segments[$segment]=array();
            $segments[$segment][$count[$segment]]=array("Type"=>$info[1], "Status"=>$info[2], "Chimeric"=>$info[3], "AssemblerID"=>$info[4]);
            $sequence="";
        }
        while (($line = fgets($handle)) !== false) {
            if (strpos($line, ">")>-1){
                $segments[$segment][$count[$segment]]['Sequence']=bzcompress($sequence,9);
                $sequence="";
                if(!strpos($line, "complete")){
                    $complete=0;
                    #echo $complete;
                }
                if (array_key_exists(intval($line[1]), $count)){
            	     $count[intval($line[1])]+=1;
            	}
            	else {
                	$count[intval($line[1])]=1;
            	}
		$segment=$line[1];
                $info=explode("|", $line);
                $segments[$segment][$count[$segment]]=array("Type"=>$info[1], "Status"=>$info[2], "Chimeric"=>$info[3], "AssemblerID"=>$info[4]);
            }
            else {
            $sequence=$sequence.$line;}
        }
        $segments[$segment][$count[$segment]]['Sequence']=bzcompress($sequence,9);

	//** check that there are 8 segments present, otherwise set complete flag to false (0) **//
        if(count($count)<8){
            $complete=0;
        }
    }
    else {
    		die("error opening the file.");
	} 
	fclose($handle);

	foreach($segments as $key=>$val0){
	#echo $key."\n";   
	#var_dump($segments[$key]);

	$query="Select * from `Sequences` WHERE ExtractID LIKE ". $extractID . " AND `Segment` = " . $key ;
        $sth=$dbh->prepare($query);
        $sth->execute();
	#if extract/segment exists in the database
        if ($row=$sth->fetchall(PDO::FETCH_ASSOC)) {
		#if there's the same number of copies of this segment in the db as I want to update:
		if (count($row)==$count[$key]){
			$i=0;
			foreach($val0 as $val){
				$id = $row[$i]["id"];
				$i+=1;
		                $query="UPDATE `vanbah01_influenza`.`Sequences` SET Type=?, Status=?, Chimeric=?, AssemblyID=?, Sequence=?, Timestamp=CURRENT_TIMESTAMP  WHERE `id`=". $id;
            			$sth=$dbh->prepare($query);
            			$sth->bindParam(1,$val["Type"]);
            			$sth->bindParam(2,$val["Status"]);
            			$sth->bindParam(3,$val["Chimeric"]);
            			$sth->bindParam(4,$val["AssemblerID"]);
            			$sth->bindParam(5,$val["Sequence"]);
            			$ro=$sth->execute();
			}

		}
		#if different number of copies of this segment in the database from what I want to update
		#just delete what's in the database and insert the new version(s)
		else {
		#delete
			foreach ($row as $key_ => $val){
			$id = $row[$key_]["id"];
			$query = "DELETE FROM `vanbah01_influenza`.`Sequences` WHERE `id` LIKE " . $id; 
                        $sth=$dbh->prepare($query);
			$sth->execute();
		}	

		#insert
			foreach($val0 as $val){
	                        $query="INSERT INTO `vanbah01_influenza`.`Sequences` (`id`, `ExtractID`, `Segment`, `Type`, `Status`, `Chimeric`, `AssemblyID`, `Sequence`, `Timestamp`) VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP);";
        	                $sth=$dbh->prepare($query);
                	        $sth->bindParam(1,$extractID);
                       		$sth->bindParam(2,$key);
	                        $sth->bindParam(3,$val["Type"]);
        	                $sth->bindParam(4,$val["Status"]);
                	        $sth->bindParam(5,$val["Chimeric"]);
                        	$sth->bindParam(6,$val["AssemblerID"]);
                        	$sth->bindParam(7,$val["Sequence"]);
                        	$ro=$sth->execute();
			}

		
		}


	}
	#if not in the database yet, insert as many versions of this segment as there's in the fasta file
	else {
                foreach($val0 as $val){
            		$query="INSERT INTO `vanbah01_influenza`.`Sequences` (`id`, `ExtractID`, `Segment`, `Type`, `Status`, `Chimeric`, `AssemblyID`, `Sequence`, `Timestamp`) VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP);";
            		$sth=$dbh->prepare($query);
            		$sth->bindParam(1,$extractID);
            		$sth->bindParam(2,$key);
            		$sth->bindParam(3,$val["Type"]);
            		$sth->bindParam(4,$val["Status"]);
            		$sth->bindParam(5,$val["Chimeric"]);
            		$sth->bindParam(6,$val["AssemblerID"]);
            		$sth->bindParam(7,$val["Sequence"]);
            		$ro=$sth->execute();
                }


	}

	/*
        $query="Select * from `Sequences` WHERE ExtractID LIKE ". $extractID . " AND `Segment` = " . $key ;
        $sth=$dbh->prepare($query);
        $sth->execute();
        if ($row=$sth->fetch(PDO::FETCH_ASSOC)) {
	    $query="UPDATE `vanbah01_influenza`.`Sequences` SET Type=?, Status=?, Chimeric=?, AssemblyID=?, Sequence=?, Timestamp=CURRENT_TIMESTAMP  WHERE ExtractID=". $extractID ." AND `Segment` = ". $key;
            $sth=$dbh->prepare($query);
	    $sth->bindParam(1,$val["Type"]);
            $sth->bindParam(2,$val["Status"]);
            $sth->bindParam(3,$val["Chimeric"]);
            $sth->bindParam(4,$val["AssemblerID"]);
            $sth->bindParam(5,$val["Sequence"]);
            #$ro=$sth->execute();
	    #echo "update segment". $key .": " .$ro."\n";
        }
        else {
            $query="INSERT INTO `vanbah01_influenza`.`Sequences` (`id`, `ExtractID`, `Segment`, `Type`, `Status`, `Chimeric`, `AssemblyID`, `Sequence`, `Timestamp`) VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP);";
	    $sth=$dbh->prepare($query);
            $sth->bindParam(1,$extractID);
            $sth->bindParam(2,$key);
            $sth->bindParam(3,$val["Type"]);
            $sth->bindParam(4,$val["Status"]);
            $sth->bindParam(5,$val["Chimeric"]);
            $sth->bindParam(6,$val["AssemblerID"]);
            $sth->bindParam(7,$val["Sequence"]);
            #$ro=$sth->execute();
            #echo "insert segment". $key .": " .$ro."\n";
        }*/
	}


	/*
	* Process Assemblies Table
	*/


    $filename=$argv[1]."_curated.report.pdf";
    $report = base64_encode(fread(fopen($filename, "r"),filesize($filename)));
    $filename=$argv[1]."_curated.variants.calls.txt";
    $variants = base64_encode(fread(fopen($filename, "rb"),filesize($filename)));
    $filename=$argv[1]."_protein.fa";
    $protein = base64_encode(fread(fopen($filename, "rb"),filesize($filename)));
    $filename=$argv[1].".anno";
    $annotations = base64_encode(fread(fopen($filename, "rb"),filesize($filename)));

    $query="Select `idAssemblies` from `Assemblies` WHERE ExtractID LIKE ". $extractID  ;
    $sth=$dbh->prepare($query);
    $ro=$sth->execute();
    if ($row=$sth->fetch()) {
	$query="UPDATE `vanbah01_influenza`.`Assemblies` SET `Fasta`='".$fasta."', `proteinFasta`='".$protein."', `annotations`='".$annotations."', `PDFReport`='".$report."', `Variants`='".$variants."', `Complete`=".$complete.", `Timestamp`=CURRENT_TIMESTAMP  WHERE ExtractID=". $extractID;
	$sth=$dbh->prepare($query);
        $ro=$sth->execute();
	echo "updating assembly:" . $ro ."\n";
    }
    else {
        $dbh->setAttribute( PDO::ATTR_ERRMODE, PDO::ERRMODE_WARNING );
        $query="INSERT INTO `vanbah01_influenza`.`Assemblies` (`idAssemblies`, `ExtractID`,`Fasta`, `proteinFasta`, `annotations`, `PDFReport`,`Variants`,`Complete`, `TimeStamp`) VALUES (NULL,".$extractID.",'".$fasta."','".$protein."','".$annotations."','".$report."','".$variants."',".$complete.",CURRENT_TIMESTAMP)";
        $sth=$dbh->prepare($query);
        $ro = $sth->execute();
	#echo $ro . "\n";
        if (!$ro) print_r($sth->errorInfo());
        echo "inserting assembly:" . $ro."\n";
    }


    if ($complete==0){
        $status="Partial";
    }
    else {
        $status="Complete";
    }

    //set subtype in extracts table according to the Segment 4 H type and segment 6 N type:
    $query="Select * from `Extracts` WHERE idExtract LIKE ". $extractID;
    $sth=$dbh->prepare($query);
    $sth->execute();
    if ($row=$sth->fetch(PDO::FETCH_ASSOC)) {
        echo "Extract exists\n";
	$subtype="";
        if (array_key_exists(4, $segments)){
		$h1=explode("N", $segments[4]["Type"]);
        	$h2=explode("H", $h1[0]);
        	$subtype=$subtype . "H" . $h2[1];
        }
	else {
	$subtype= $subtype."na";
	}
	 if (array_key_exists(6,$segments)){
		$n1=explode("N", $segments[6]["Type"]);
        	$subtype= $subtype."N" . $n1[1];
        }
	else {
	$subtype= $subtype . "na";
	}
	#$subtype="H".$h."N".$n;
        $query="UPDATE `vanbah01_influenza`.`Extracts` SET `Status` = ?, `Subtype`=?, Modified=CURRENT_TIMESTAMP WHERE `idExtract` = ". $extractID;
        $sth=$dbh->prepare($query);
        $sth->bindParam(1,$status);
        $sth->bindParam(2,$subtype);
        $ro=$sth->execute();
        echo "update extract table Subtype with type ".$subtype .": " . $ro . "\n";
    }



}

$dbh=null;

