Sub ºê1()
'
' ºê1 ºê
'

'
    Dim C As String
    Dim n As Integer
    For i = 0 To 12
        n = i * 5 + 20
        Windows("concrete_t.xls").Activate
        Sheets("C" & CStr(n)).Select
        Rows("1:1").Select
        Selection.Insert Shift:=xlDown, CopyOrigin:=xlFormatFromLeftOrAbove
        Selection.Insert Shift:=xlDown, CopyOrigin:=xlFormatFromLeftOrAbove
        
        Range("A1:B1").Select
        With Selection
            .HorizontalAlignment = xlCenter
            .VerticalAlignment = xlCenter
            .WrapText = False
            .Orientation = 0
            .AddIndent = False
            .IndentLevel = 0
            .ShrinkToFit = False
            .ReadingOrder = xlContext
            .MergeCells = False
        End With
        Selection.Merge
        Range("A2").Select
        ActiveCell.FormulaR1C1 = "Ó¦±ä"
        Range("B2").Select
        ActiveCell.FormulaR1C1 = "Ó¦Á¦"
        Range("A1:B1").Select
        ActiveCell.FormulaR1C1 = "ÊÜÀ­Mpa"
        
        Range("C1:D1").Select
        With Selection
            .HorizontalAlignment = xlCenter
            .VerticalAlignment = xlCenter
            .WrapText = False
            .Orientation = 0
            .AddIndent = False
            .IndentLevel = 0
            .ShrinkToFit = False
            .ReadingOrder = xlContext
            .MergeCells = False
        End With
        Selection.Merge
        Range("C2").Select
        ActiveCell.FormulaR1C1 = "Ó¦±ä"
        Range("D2").Select
        ActiveCell.FormulaR1C1 = "Ó¦Á¦"
        Range("C1:D1").Select
        ActiveCell.FormulaR1C1 = "ÊÜÑ¹Mpa"
        
        Windows("concrete_c.xls").Activate
        Sheets("C" & CStr(n)).Select
        ActiveWindow.SmallScroll Down:=-21
        Range("A1:B1").Select
        Range(Selection, Selection.End(xlDown)).Select
        Selection.Copy
        Windows("concrete_t.xls").Activate
        Sheets("C" & CStr(n)).Select
        Range("C3").Select
        ActiveSheet.Paste
        ActiveWindow.SmallScroll Down:=57
    
        Next i
    End
    
End Sub