CREATE TABLE Raft (
   id int NOT NULL AUTO_INCREMENT,
   raftId char(2) NOT NULL,
   PRIMARY KEY (id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE CCD (
   id int NOT NULL AUTO_INCREMENT,
   raftId int,
   raftSlot char(2),
   gainMedian float,
   fullWellMean float,   
   fullWellMedian float,   
   fullWellStdev float,
   maxDeviation float,
   ctiSerialMean float,
   ctiParallelMean float,
   darkCurrent95Mean float,
   numBrightPixels float,
   PRIMARY KEY (id),
   CONSTRAINT fk20 FOREIGN KEY (raftId)
   REFERENCES Raft (id), INDEX fk20 (raftId)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE Segment (
   id int NOT NULL AUTO_INCREMENT,
   channelId char(2) NOT NULL,
   ccdId int NOT NULL,
   gain float,
   readNoise float,
   linefit_Slope float,
   linefit_Intercept float,
   maxDeviation float,
   fullWell float,
   ctiSerial float,
   ctiParallel float,
   darkCurrent95 float,
   numBrightPixels int,
   PRIMARY KEY (id),
   CONSTRAINT fk10 FOREIGN KEY (ccdId)
   REFERENCES CCD (id), INDEX fk10 (ccdId)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE CCD_VendorIds (
   id int NOT NULL AUTO_INCREMENT,
   ccdId int NOT NULL,
   vendor varchar(30) NOT NULL,
   vendorId varchar(30) NOT NULL,
   PRIMARY KEY (id),
   CONSTRAINT fk30 FOREIGN KEY (ccdId)
   REFERENCES CCD (id), INDEX fk30 (ccdId),
   UNIQUE (vendor, vendorId)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
