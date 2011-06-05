
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->



<p> <strong>Abstract:</strong>  Constructing phylogenetic trees from measures of dissimilarity (distances) between species is a popular approach for inferring phylogenies, as it is fast and provides a useful snapshot of the data. The most popular method for such reconstruction is Neighbour Joining (N.J). Many improvements have been made to this algorithm since its release, for instance Bio-NJ, which have sped NJ up and have improved its accuracy. However in practice distance matrices with missing values (due to noise or lack of confidence in readings) are encountered. For such cases modifications of the above methods have been made, such as NJ*, which given an incomplete distance matrix will return a phylogenetic tree from that distance matrix. The project aims to extends the APE R package (which currently contains implementations of NJ and Bio_NJ) to also include such methods for phylogenetic reconstruction. We will also add a non-agglomerative approach (the triangles method) to APE and later extend this to handle reconstruction from incomplete distances.
</p>

<p> <strong>Project Plan</strong> <br>

<bf>

</p>




<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>



</body>

</html>
