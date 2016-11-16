<html>
    <head>
        <title>LINbase</title>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.2/jquery.min.js"></script>
        <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">
        <link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">
        <script type = "text/javascript" src="<?php echo base_url().'CodeIgniter/JS/functions.js'; ?>"></script>
        <link rel="stylesheet" type="text/css" href="//cdn.datatables.net/1.10.12/css/jquery.dataTables.css">
        <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/grids-responsive-min.css">
        <link rel="stylesheet" type="text/css" href="<?php echo base_url().'CodeIgniter/CSS/generalStyle.css'; ?>">
    </head>
<body>
<div class="container" style="width:1200px;">
    <div id="center">
        <div>
            <h1 style="text-align: center">LINbase: The Life Identification Numbers Platform</h1>
            <label><h3 style="font-weight: 100;color: rgba(0, 0, 0, 0.27);text-align: right">Where you share everything about microorganisms.</h3></label>
        </div>
        <div style="text-align: right; margin-right: 25px;">
            <a href="#" id="logout">Logout</a>
        </div>
        <br><br>
        <?php
        $servername = "128.173.74.68";
        $username = "linproject";
        $password = "tcejorpnil";
        $db_name = "LINdb_Psy";
        $conn = new mysqli($servername, $username, $password, $db_name);
        if ($conn->connect_error) {
            die("Connection failed: " . $conn->connect_error);
        }
        session_start();
        $LIN = $_SESSION["LIN"];
        $Description_Item_Names = $_SESSION["Description_Item_Name"];
        $Description_Item_Names_chosen = [];
        $description_option = []; // Description_Item_ID chosen
        $description_text = []; // DescriptionValue entered
        //$row_check = []; // Index of checked row
        $col_check = []; // Indices of checked columns, consecutive and starts from 0

        for($i=0;$i<count($_POST["description_input"]);$i++){
            if($_POST["description_input"][$i] != 0 && $_POST["description_input_text"] != ''){
                $description_option[] = $_POST["description_input"][$i];
                $description_text[] = $_POST["description_input_text"][$i];
                $sql = "SELECT Description_Item_Name FROM Description_Items WHERE Description_Item_ID={$_POST["description_input"][$i]}";
                $result = $conn->query($sql);
                $Description_Item_Names_chosen[] = $result->fetch_assoc()["Description_Item_Name"];
            }
        }

        if(!empty($_POST["row_checkbox"])){
            $row_check = array_shift($_POST["row_checkbox"]);
        }

        if(!empty($_POST["col_checkbox"])){
            foreach($_POST["col_checkbox"] as $each_col_checkbox){
                $col_check[] = $each_col_checkbox;
            }
        }

        $LIN_to_describe_whole = $LIN[$row_check];
        $parsed_LIN_to_describe_whole = explode(",", $LIN_to_describe_whole);
        $parsed_LIN_to_describe_part = array_slice($parsed_LIN_to_describe_whole, 0, count($col_check));

        $LIN_to_describe_part = implode(",",$parsed_LIN_to_describe_part);
        $LIN_to_describe_part = "'{$LIN_to_describe_part}'";
        $Load_Description = [];
        for($i=0;$i<count($description_option);$i++){
            $sql = "
            INSERT INTO Description (Part_LIN, Description_Item_ID, DescriptionValue) VALUES ($LIN_to_describe_part,
            {$description_option[$i]}, '{$description_text[$i]}')
            ";
            $Load_Description[] = $sql;
        }

        $description_items = implode(",",$Description_Item_Names_chosen);
        $description_items_to_enter = implode("','",$description_text);
        $description_items_to_enter = "'{$description_items_to_enter}'";
        $Load_LIN_to_Description = "
        INSERT INTO LIN_to_Description (Part_LIN, {$description_items}) VALUES ({$LIN_to_describe_part}, {$description_items_to_enter})
        ";
        $_SESSION["Load_Description"] = $Load_Description;
        $_SESSION["Load_LIN_to_Description"] = $Load_LIN_to_Description;
        ?>
        <form action="submit_description.php" method="post" class="pure-form pure-form-aligned">
            <fieldset>
                <legend>Confirm your description to add</legend>
                <table>
                    <thead>
                    <tr>
                        <th>LINs to describe</th>
                        <?php
                        foreach ($Description_Item_Names_chosen as $each_chosen_description){
                            echo "<th>{$each_chosen_description}</th>";
                        }
                        ?>
                    </tr>
                    </thead>
                    <tbody>
                    <tr>
                        <?php
                        echo "<td>{$LIN_to_describe_part}</td>";
                        foreach ($description_text as $each_description_text){
                            echo "<td>{$each_description_text}</td>";
                        }
                        ?>
                    </tr>
                    </tbody>
                </table>
            </fieldset>
            <br><br>
            <div style="text-align: right">
                <button class="pure-button pure-button-primary" type="submit">Confirm</button>
            </div>
        </form>
    </div>
</div>
    <footer class="footer">
			<div class="container">
				<footer>
					The Life Identification Numbers Platform is co-developed by
					the Department of Plant Pathology, Physiology and Weed Science and
					the Department of Computer Science of Virginia Polytechnic Institute and State University.
					For any questions or problems, please contact Long Tian by longtian@vt.edu .
				</footer>
			</div>
		</footer>
</body>
</html>